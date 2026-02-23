import os
import pathlib
import tkinter as tk
from tkinter import filedialog

def delete_temp_files(root_dir):
    root_path = pathlib.Path(root_dir)
    deleted = []

    print(f"\nSearching in: {root_path}")
    print("Looking for *.bin and temp_wh.dat files...\n")

    # First collect all matching files (so we can report how many)
    matches = []
    for path in root_path.rglob("*"):
        if path.is_file() and (path.suffix == ".bin" or path.name == "temp_wh.dat"):
            matches.append(path)

    if not matches:
        print("No matching temporary files found. Nothing to delete.")
        return

    print(f"Found {len(matches)} matching files. Starting deletion...\n")

    # Delete with per-file messages
    for i, path in enumerate(matches, start=1):
        print(f"[{i}/{len(matches)}] Deleting: {path}")
        path.unlink()
        deleted.append(str(path))

    print("\nDeleted files:")
    for f in deleted:
        print(f)
    print(f"\nTotal deleted: {len(deleted)}")

if __name__ == "__main__":
    # Ask user to pick a root folder
    tk_root = tk.Tk()
    tk_root.withdraw()  # hide main window
    folder = filedialog.askdirectory(title="Select root folder")
    tk_root.destroy()

    if folder:
        print(f"Selected folder: {folder}")
        delete_temp_files(folder)
    else:
        print("No folder selected, exiting.")
