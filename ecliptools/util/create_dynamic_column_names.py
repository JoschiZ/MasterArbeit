#!/usr/bin/python
from pathlib import Path


def create_dynamic_column_names(file_path: Path, delimiter="\t"):
    # The max column count a line in the file could have
    largest_column_count = 0

    # Loop the data lines
    with open(file_path, 'r') as temp_f:
        # Read the lines
        lines = temp_f.readlines()

        for line in lines:
            # Count the column count for the current line
            column_count = len(line.split(delimiter)) + 1

            # Set the new most column count
            largest_column_count = column_count if largest_column_count < column_count else largest_column_count

    # Generate column names (will be 0, 1, 2, ..., largest_column_count - 1)
    column_names = [column for column in range(0, largest_column_count)]
    return column_names

