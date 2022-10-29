# fastz

Can I buy this molecule? Returns results in about 500 ns and consumes about 1 GB of memory.

```py

from fastz import purchasable
purchasable('CCCO')

```

If it returns `True` - it may be purchasable. If it returns `False` - it is not purchasable. Reference is ZINC20 October batch. Add `canonicalize=True` if your SMILES are not canonicalized.