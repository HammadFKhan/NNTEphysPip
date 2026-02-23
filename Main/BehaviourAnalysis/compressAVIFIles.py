from pathlib import Path
from moviepy.editor import VideoFileClip
import tkinter as tk
from tkinter import filedialog

def compress_avi_to_mp4(input_path, output_path,
                        resolution=None,
                        bitrate="3000k"):
    """
    Compress a single AVI to MP4 (H.264) using MoviePy.

    input_path: Path to .avi file
    output_path: Path to .mp4 file
    resolution: (width, height) or None to keep original
    bitrate: target video bitrate string, e.g. "2000k"
    small -> smaller files, lower quality; 
    higher → larger files, better quality.
    """
    print(f"\nLoading: {input_path}")
    clip = VideoFileClip(str(input_path))

    if resolution is not None:
        clip = clip.resize(resolution)

    clip.write_videofile(
        str(output_path),
        codec="libx264",      # H.264
        audio_codec="aac",    # reasonable default
        bitrate=bitrate
    )

def batch_compress_folder(root_dir,
                          resolution=None,
                          bitrate="3000k"):
    root_path = Path(root_dir)
    avi_files = list(root_path.rglob("*.avi"))

    if not avi_files:
        print("No .avi files found.")
        return

    print(f"Found {len(avi_files)} .avi files under {root_path}\n")

    for i, avi_path in enumerate(avi_files, start=1):
        mp4_path = avi_path.with_suffix(".mp4")  # same name, .mp4 extension
        print(f"[{i}/{len(avi_files)}] Compressing {avi_path.name} -> {mp4_path.name}")
        compress_avi_to_mp4(avi_path, mp4_path,
                            resolution=resolution,
                            bitrate=bitrate)

if __name__ == "__main__":
    # Ask user to pick a root folder
    tk_root = tk.Tk()
    tk_root.withdraw()
    folder = filedialog.askdirectory(title="Select root folder containing .avi files")
    tk_root.destroy()

    if folder:
        print(f"Selected folder: {folder}")
        # Adjust resolution/bitrate as needed, or set resolution=None to keep original size
        batch_compress_folder(folder,  
                              bitrate="3000k") 
    else:
        print("No folder selected, exiting.")
