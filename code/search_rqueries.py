import os
import re
import json
import sys

def binary_search_json(file_path, target_key:int, chunk_size=1024):
    """
    Searches for a key in a massive JSON file using dichotomy.
    Assumes the file is a single flat dictionary with sorted keys.

    The chunk size must be large enough to hold any key.
    """
    file_size = os.path.getsize(file_path)
    low = 0
    high = file_size

    # Regex to find a key boundary:
    key_pattern = re.compile(rb'"(\d+)"\s*:')

    with open(file_path, 'rb') as f:
        while low <= high:
            mid = (low + high) // 2

            f.seek(mid)
            chunk = f.read(chunk_size)

            if not chunk:
                # Hit EOF without finding anything; target must be before mid
                high = mid - 1
                continue

            matches = list(key_pattern.finditer(chunk))

            valid_match = None
            for m in matches:
                try:
                    current = int(m.group(1))
                    valid_match = (m, current)
                    break
                except ValueError:
                    continue

            if not valid_match:
                # No key found after 'mid' in this chunk. This usually means
                # we landed inside a massive value or near the end of the file.
                # In either case, the target must be before 'mid'.
                high = mid - 1
                continue

            match_obj, current = valid_match
            key_absolute_offset = mid + match_obj.start()

            if current == target_key:
                # KEY FOUND!
                # Move the cursor to right after the colon to parse the value
                value_start_offset = mid + match_obj.end()
                f.seek(value_start_offset)

                decoder = json.JSONDecoder()
                buffer = ""
                read_size = 4096

                while True:
                    # Read the file in chunks until we have a JSON value
                    data = f.read(read_size)
                    if not data:
                        # Reached EOF unexpectedly
                        return None

                    buffer += data.decode('utf-8', errors='replace')

                    try:
                        # raw_decode extracts exactly one complete JSON value (string, dict, list, etc.)
                        val, _ = decoder.raw_decode(buffer.lstrip())
                        return val
                    except json.JSONDecodeError:
                        # The buffer doesn't contain a complete JSON value yet; keep reading
                        pass

            elif current < target_key:
                # The key we found is smaller than the target.
                # The target must exist strictly after this key.
                low = key_absolute_offset + 1
            else:
                # The first key we found AFTER 'mid' is greater than the target.
                # Because it was the *first* one found, the target must be BEFORE 'mid'.
                high = mid - 1

    return None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Usage: search_rqueries.py <file> <key1> <key2> [...]")

    file_name, *keys = sys.argv[1:]

    for k in keys:
        result = binary_search_json(file_name, int(k))
        if result is not None:
            print(f"{k}: {result}")
        else:
            print(f"{k}: NOT FOUND")
