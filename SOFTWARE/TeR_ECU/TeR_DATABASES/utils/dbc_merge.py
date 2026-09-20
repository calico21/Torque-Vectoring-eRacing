#!/usr/bin/env python3
#Ozuba 2024

import argparse
import cantools

def merge_dbc(dbc_files, output_file):
    """
    Merge multiple CAN DBC databases into one and save to a file.

    Args:
        dbc_files (list of str): List of paths to DBC files.
        output_file (str): Path to save the merged DBC file.
    """
    merged_db = cantools.database.Database()
    for dbc_file in dbc_files:
        db = cantools.database.load_file(dbc_file)
        for message in db.messages:
            # Check for ID conflicts
            if any(msg.frame_id == message.frame_id for msg in merged_db.messages):
                print(f"Conflict detected for message ID {message.frame_id:#x}. Skipping {message.name}.")
                continue
            merged_db.messages.append(message)

        # Optionally merge nodes or other definitions
        for node in db.nodes:
            if node not in merged_db.nodes:
                merged_db.nodes.append(node)

    # Save the merged database
    with open(output_file, 'w') as f:
        f.write(merged_db.as_dbc_string())

    print(f"Merged DBC saved to '{output_file}'.")

def main():
    parser = argparse.ArgumentParser(
        description="Merge multiple DBC files into a single DBC file."
    )
    parser.add_argument(
        "dbc_files",
        nargs='+',
        help="Paths to the DBC files to merge.",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to save the merged DBC file.",
    )

    args = parser.parse_args()

    try:
        merge_dbc(args.dbc_files, args.output)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
