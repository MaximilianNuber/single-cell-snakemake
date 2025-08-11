from __future__ import annotations
import inspect

class ContractError(TypeError):
    pass

def _required_positional_params(fn) -> int:
    """
    Count required positional-or-keyword parameters on a *bound* method.
    ('self' is already bound and not counted.)
    """
    sig = inspect.signature(fn)
    count = 0
    for p in sig.parameters.values():
        if p.kind in (inspect.Parameter.POSITIONAL_ONLY,
                      inspect.Parameter.POSITIONAL_OR_KEYWORD):
            if p.default is inspect._empty:
                count += 1
    return count

def ensure_methods(obj, *, require: dict[str, int]) -> None:
    """
    Ensure `obj` has each method in `require` and that each method has at least
    the given number of required positional parameters.

    Example: ensure_methods(emb, require={"compute": 3})
    """
    cls = obj.__class__.__name__
    for name, min_required in require.items():
        fn = getattr(obj, name, None)
        if fn is None or not callable(fn):
            raise ContractError(f"{cls} is missing a callable method '{name}()'.")
        req = _required_positional_params(fn)
        if req < min_required:
            raise ContractError(
                f"{cls}.{name}() requires {req} positional args, "
                f"but at least {min_required} are expected."
            )