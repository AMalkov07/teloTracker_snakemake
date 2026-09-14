# Reads the `cut99` scheme gets wrong

18 misread reads over 2 day-0 populations; 18 wrong Y' copies.

## What is mislabelled as what

| true element | matched element | true group | assigned group | copies | reads | chr ends | samples |
|---|---|---|---|---|---|---|---|
| chr6L-1 | chr14L-5 | `T99-1` | `T99-2` | 4 | 4 | chr6L (4) | 1 |
| chr16R-1 | chr14L-3 | `T99-6` | `T99-2` | 2 | 2 | chr16R (2) | 1 |
| chr14R-1 | chr14L-1 | `T99-10` | `T99-8` | 2 | 2 | chr14R (2) | 1 |
| chr5R-1 | chr13L-1 | `T99-9` | `T99-2` | 1 | 1 | chr5R (1) | 1 |
| chr6L-1 | chr14L-3 | `T99-1` | `T99-2` | 1 | 1 | chr6L (1) | 1 |
| chr10L-1 | chr14R-1 | `T99-7` | `T99-10` | 1 | 1 | chr10L (1) | 1 |
| chr2L-1 | chr8R-1 | `T99-1` | `T99-3` | 1 | 1 | chr2L (1) | 1 |
| chr5R-1 | chr10L-1 | `T99-9` | `T99-7` | 1 | 1 | chr5R (1) | 1 |
| chr5R-1 | chr14R-1 | `T99-9` | `T99-10` | 1 | 1 | chr5R (1) | 1 |
| chr6L-1 | chr16L-1 | `T99-1` | `T99-8` | 1 | 1 | chr6L (1) | 1 |
| chr6L-1 | chr14R-1 | `T99-1` | `T99-10` | 1 | 1 | chr6L (1) | 1 |
| chr6L-1 | chr8R-1 | `T99-1` | `T99-3` | 1 | 1 | chr6L (1) | 1 |
| chr8R-1 | chr2L-1 | `T99-3` | `T99-1` | 1 | 1 | chr8R (1) | 1 |

## The reads

| sample | end | read | copies | wrong at | expected elements | matched elements |
|---|---|---|---|---|---|---|
| 7302_day0_with_selection | chr10L | `30110432-fc36-40ad-a8e5-b0520fcb16a2` | 1 | 1 | chr10L-1 | chr14R-1 |
| 7302_day0_with_selection | chr14R | `4bf91fa6-c82f-4385-a97f-e28e1ae4a55f` | 1 | 1 | chr14R-1 | chr14L-1 |
| 7302_day0_with_selection | chr14R | `5c5de926-93b8-4eb1-bfb7-854b242ef9bd` | 1 | 1 | chr14R-1 | chr14L-1 |
| 7172_day0_with_selection | chr16R | `47be2d6c-e32d-4581-b8a4-0097d68fbbd0` | 1 | 1 | chr16R-1 | chr14L-3 |
| 7172_day0_with_selection | chr16R | `856225af-6cb4-45e3-b958-0178ec1f0eb5` | 1 | 1 | chr16R-1 | chr14L-3 |
| 7302_day0_with_selection | chr2L | `d8f7849a-17d8-4f31-992f-64fd4334c2aa` | 1 | 1 | chr2L-1 | chr8R-1 |
| 7172_day0_with_selection | chr5R | `40b0a0d2-0382-4db6-8f84-dd4a0caa8813` | 1 | 1 | chr5R-1 | chr13L-1 |
| 7302_day0_with_selection | chr5R | `b38bcd14-a668-472b-9de2-a01b5c1036c2` | 1 | 1 | chr5R-1 | chr10L-1 |
| 7302_day0_with_selection | chr5R | `214d0403-7b3d-4f40-89e1-61984253d8ba` | 1 | 1 | chr5R-1 | chr14R-1 |
| 7172_day0_with_selection | chr6L | `5c6570e5-11da-41bd-b830-a6513179c7db` | 1 | 1 | chr6L-1 | chr14L-3 |
| 7302_day0_with_selection | chr6L | `a3078675-718a-42df-a6da-b0e8306e4ba2` | 1 | 1 | chr6L-1 | chr16L-1 |
| 7302_day0_with_selection | chr6L | `3a398ad6-3c1a-41a1-ae47-a3e357022d1d` | 1 | 1 | chr6L-1 | chr14R-1 |
| 7302_day0_with_selection | chr6L | `6d06047e-4254-40f4-ae6c-e698a015fa1c` | 1 | 1 | chr6L-1 | chr14L-5 |
| 7302_day0_with_selection | chr6L | `a10ea8e7-b6e6-4544-9629-458d3169ef38` | 1 | 1 | chr6L-1 | chr14L-5 |
| 7302_day0_with_selection | chr6L | `a22f7319-a8b1-4e2d-95f6-6e0df7d1cefb` | 1 | 1 | chr6L-1 | chr14L-5 |
| 7302_day0_with_selection | chr6L | `fdc4737e-963b-4563-907f-0d7174c68d67` | 1 | 1 | chr6L-1 | chr8R-1 |
| 7302_day0_with_selection | chr6L | `a7d76bd0-ec4f-4ba4-88d7-388bd3cde346` | 1 | 1 | chr6L-1 | chr14L-5 |
| 7302_day0_with_selection | chr8R | `1b96eb7b-781f-43e8-9b6b-9d80dc4591f3` | 1 | 1 | chr8R-1 | chr2L-1 |
