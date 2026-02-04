import sys
from pathlib import Path

def add_src_to_path(package_dirname: str = "rayleigh_gans") -> Path:
    """
    Add repo_root/src to sys.path by searching upward for src/<package_dirname>.
    Returns the detected repo_root.
    """
    p = Path.cwd().resolve()
    for parent in [p] + list(p.parents):
        if (parent / "src" / package_dirname).exists():
            src_path = str(parent / "src")
            if src_path not in sys.path:
                sys.path.insert(0, src_path)
            return parent
    raise RuntimeError(
        f"Could not find src/{package_dirname} from {p} or its parents."
    )