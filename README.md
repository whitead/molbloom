# fastz

Can I buy this molecule? Results in 0.1ms or your money back.

```py

from fastz import purchasable
purchasable('CCCO')

```

If it returns `True` - it may be purchasable. If it returns `False` - it is not purchasable. Reference is ZINC20 October batch. Add `canonicalize=True` if your SMILES are not canonicalized.