#!/usr/bin/env python3
import csv

# ===== 只改这里 =====
output_file = "channel_map.csv"

# board_id -> detector_id
board_det = {
    0: 1,
    2: 2,
    3: 3,
}

pedestal = 0
noise_sigma = 0
gain = 1
polarity = -1
status = 0
# ===================


def apv_channel_to_strip(ch):
    return 127 - ch // 16 - 8 * ((ch % 16) // 4) - 32 * (ch % 4)


header = [
    "board_id",
    "chip_id",
    "channel_id",
    "detector_id",
    "plane_type",
    "strip_id",
    "polarity",
]

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)

    for board_id, detector_id in board_det.items():
        for chip_id in range(4):
            # chip 0,1 -> plane 0
            # chip 2,3 -> plane 1
            plane_type = chip_id // 2

            # chip 0,2 -> strip 0-127
            # chip 1,3 -> strip 128-255
            strip_offset = 128 * (chip_id % 2)

            for channel_id in range(128):
                strip_id = strip_offset + apv_channel_to_strip(channel_id)

                writer.writerow([
                    board_id,
                    chip_id,
                    channel_id,
                    detector_id,
                    plane_type,
                    strip_id,
                    polarity,
                ])

print("Wrote", output_file)