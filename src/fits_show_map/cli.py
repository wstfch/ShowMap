from __future__ import annotations

import argparse

from .core import show_map


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fits-show-map",
        description="Display a FITS image with WCS axes and optional PNG export.",
    )
    p.add_argument("image", help="Path to a FITS file")
    p.add_argument("-c", "--cmap", default="viridis", help="Matplotlib colormap name")
    p.add_argument(
        "-s",
        "--savefig",
        action="store_true",
        help="If set, save a PNG next to the FITS file (same basename).",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    # Keep compatibility with the original implementation, which expects 'True'/'False'
    show_map(image=args.image, cmap=args.cmap, savefig="True" if args.savefig else "False")


if __name__ == "__main__":
    main()
