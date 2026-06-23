# Channel map format

CSV columns:

```text
board_id,chip_id,channel_id,detector_id,plane_type,strip_id,pedestal,noise_sigma,gain,polarity,status
```

The first seven columns are required. Calibration and status columns are optional.
`plane_type` uses the BeamAnalysis detector convention, normally `0` for X and
`1` for Y.
