#!/usr/bin/env python3
"""
name_to_graph.py - Convert chemical names to SMILES notation using Leruli API

This script reads chemical names from a file, sends each name to the Leruli API,
and saves the resulting SMILES strings to a new file.
"""

import asyncio
import argparse
import httpx
from pathlib import Path
from tqdm import tqdm
from tenacity import (
    retry,
    stop_after_attempt,
    wait_exponential,
    retry_if_exception_type,
)
from google import genai
from google.genai import types
from molbloom import canon

# Initialize client
client = genai.Client(vertexai=True, project="paperqa-385600", location="us-west1")


async def get_smiles(reagent):
    """Async version to get target information with caching"""

    # Otherwise, call the API
    model = "gemini-2.0-flash-001"
    contents = [
        types.Content(
            role="user",
            parts=[
                types.Part.from_text(
                    text=(
                        f"Provide SMILES for this organic chemistry reagent {reagent}."
                        " If it refers to a class, list one per line in the response."
                        " Only provide SMILES, no other information."
                    )
                ),
            ],
        )
    ]

    tools = [types.Tool(google_search=types.GoogleSearch())]

    generate_content_config = types.GenerateContentConfig(
        temperature=0.5,
        top_k=4,
        max_output_tokens=8192,
        tools=tools,
        response_mime_type="text/plain",
    )

    response = await client.aio.models.generate_content(
        model=model,
        contents=contents,
        config=generate_content_config,
    )
    return [line for line in response.text.split("\n") if line]


async def fetch_graph(client, name, semaphore):
    """
    Fetch molecular graph (SMILES) for a given chemical name.
    Uses semaphore to control concurrency.
    """
    async with semaphore:
        try:
            return name.strip(), await _fetch_with_retry(client, name.strip())
        except Exception as e:
            print(f"Error processing '{name.strip()}': {str(e)}")
            return name.strip(), None


@retry(
    retry=retry_if_exception_type((httpx.HTTPError, httpx.TimeoutException)),
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=10),
    reraise=True,
)
async def _fetch_with_retry(client, name):
    """
    Internal function to perform API requests with retries.
    Will retry up to 3 times with exponential backoff for HTTP errors.
    """
    response = await client.post(
        "https://api.leruli.com/v22_1/name-to-graph", json={"name": name}
    )
    response.raise_for_status()
    return response.json()


def extract_smiles(result):
    """
    Extract SMILES string from API response.
    The Leruli API returns a JSON object with a 'graph' field containing the SMILES notation.
    Example response: {"graph": "CC(=O)OC1=CC=CC=C1C(=O)O", "reference": "wikidata"}
    """
    if isinstance(result, dict) and "graph" in result:
        return result["graph"]
    return None


async def process_file(input_file, output_file, max_concurrency=5):
    """
    Process each line of the input file concurrently up to max_concurrency,
    and write results to the output file.
    """
    # Read all lines from input file
    with open(input_file, "r") as f:
        lines = [line for line in f if line.strip()]

    total_lines = len(lines)
    if total_lines == 0:
        print("Warning: Input file is empty.")
        return

    # Create semaphore to control concurrency
    semaphore = asyncio.Semaphore(max_concurrency)

    # Set up progress bar
    progress = tqdm(total=total_lines, desc="Processing chemical names")

    # Create output file
    with open(output_file, "w") as out_f:
        # Use a reasonable timeout for API requests
        async with httpx.AsyncClient(timeout=30.0) as client:
            tasks = []
            for line in lines:
                task = asyncio.create_task(fetch_graph(client, line, semaphore))
                task.add_done_callback(lambda _: progress.update(1))
                tasks.append(task)

            # Wait for all tasks to complete and write results
            results = await asyncio.gather(*tasks)

            success_count = 0
            error_count = 0

            for name, result in results:
                if result:
                    smiles = extract_smiles(result)
                    if smiles:
                        out_f.write(f"{canon(smiles)}\n")
                        success_count += 1
                    else:
                        print(
                            f"Warning: Could not extract SMILES for '{name}'. Response does not contain 'graph' field: {result}"
                        )
                        maybe_smiles = await get_smiles(name)
                        for s in maybe_smiles:
                            try:
                                s = canon(s)
                            except Exception:
                                print(
                                    f"Warning: Could not canonicalize SMILES for: {s}"
                                )
                                continue
                            if s:
                                out_f.write(f"{s}\n")
                            success_count += 1
                        if not maybe_smiles:
                            error_count += 1
                else:
                    error_count += 1

    progress.close()

    print(f"Processing complete: {success_count} successful, {error_count} failed.")


def main():
    """Main entry point"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Convert chemical names to SMILES notation using Leruli API"
    )
    parser.add_argument("input_file", help="Input file with one chemical name per line")
    parser.add_argument(
        "-c",
        "--concurrency",
        type=int,
        default=5,
        help="Maximum number of concurrent requests (default: 5)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="smiles.txt",
        help="Output file for SMILES strings (default: smiles.txt)",
    )

    args = parser.parse_args()

    # Verify input file exists
    if not Path(args.input_file).exists():
        print(f"Error: Input file '{args.input_file}' not found.")
        return 1

    # Run the async process
    asyncio.run(process_file(args.input_file, args.output, args.concurrency))

    print(f"Results saved to {args.output}")
    return 0


if __name__ == "__main__":
    exit(main())
