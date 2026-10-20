# SOXS UV-VIS Wavelength Calibration Analysis

This README documents the key results of adjusting and measuring the SOXS UVVIS wavelength calibration accuracy. 

* See here for the analysis and subsequent method used to [Tightening the SOXS UVVIS wavelength solution.md](Tightening%20the%20SOXS%20UVVIS%20wavelength%20solution.md).
* See here for the method used to [quantify the wavelength calibration accuracy](Measuring%20the%20UVVIS%20Wavelength%20Calibration%20Accuracy.md) (Geza's modified notebook). 

## Method to Measure Wavelength Calibration Accuracy

1. **Measure line shifts.** For each line, find the flux minimum and subtract the reference minimum. For **stellar** lines, use the star's rest-frame spectrum. **Telluric** lines use the telluric model. Each line is assigned to its order (u/g/r/i).
2. **Barycentric correction.** locations of stellar lines are corrected for Earth's velocity and rotation. Telluric lines form in Earth's atmosphere, so they stay a fixed reference.
3. **Fit star velocity.** In the r and i orders, take the median telluric offset per exposure, then compute `vStar` from each stellar line: shift minus telluric offset, converted to velocity, plus the barycentric correction. Take the 3-sigma MAD-clipped median.
4. **Compute residuals.** Subtract the expected star-velocity shift from each line's shift. Telluric lines expect 0, so their residual equals their shift. The residual is the wavelength-solution error.
5. **Per-order correction.** Compute clipped stats per order. The **u and g orders anchor on stellar lines** because they have no telluric lines. The **r and i orders anchor on telluric lines**, because their stellar lines fit the velocity and their residuals are about 0 by construction. The correction is minus the median shift.


Here are the before-and-after results of the skyline wavelength-calibration enhancement. The ‘before’ results come from data reduced with soxpipe v0.17.4, and the ‘after’ results from v0.19.0.

For each reduction, we reduce the same 11 standard stars. First, we reduce them with the daily arcs (normal operations mode), then with the attached arcs. 

## BEFORE (v0.17.4) with Daily Arcs 

<img width="1584" height="1084" alt="image" src="https://github.com/user-attachments/assets/b9698feb-6c71-4c1d-84d5-eed62a0956a9" />

RMS = √(median residual² + σ²), pooling telluric shifts and velocity-corrected stellar residuals.

| Band      |   N |   Median [Å] |   σ [Å] |   RMS [Å] |
|:----------|----:|-------------:|--------:|----------:|
| u         |  20 |       +0.007 |   0.363 |     0.363 |
| g         |  25 |       -0.209 |   0.596 |     0.632 |
| r         |  23 |       +0.277 |   0.389 |     0.478 |
| i         |  33 |       +0.186 |   0.600 |     0.628 |
| All bands | 101 |       +0.072 |   0.592 |     0.596 |

**Telluric lines in i-band**

<img width="1548" height="441" alt="image" src="https://github.com/user-attachments/assets/e441858a-4547-4ebd-9aac-fb7ec627fc94" />

**H-beta stellar line in g-band**

<img width="1542" height="745" alt="image" src="https://github.com/user-attachments/assets/fab0897d-fd65-460b-968a-7d7559aad760" />


## BEFORE (v0.17.4) with Attached Arcs 

<img width="1585" height="1081" alt="image" src="https://github.com/user-attachments/assets/8e8f1d6c-539a-4e7f-9a08-2ae858d00536" />

RMS = √(median residual² + σ²), pooling telluric shifts and velocity-corrected stellar residuals.

| Band      |   N |   Median [Å] |   σ [Å] |   RMS [Å] |
|:----------|----:|-------------:|--------:|----------:|
| u         |  21 |       +0.201 |   0.369 |     0.420 |
| g         |  25 |       +0.012 |   0.678 |     0.678 |
| r         |  23 |       +0.021 |   0.203 |     0.204 |
| i         |  33 |       -0.077 |   0.223 |     0.236 |
| All bands | 102 |       +0.007 |   0.418 |     0.418 |

**Telluric lines in i-band**

<img width="1534" height="749" alt="image" src="https://github.com/user-attachments/assets/dfab0436-c115-47a2-8158-4d2e7b6563ae" />

**H-beta stellar line in g-band**

<img width="1543" height="447" alt="image" src="https://github.com/user-attachments/assets/2da77ddd-4ecb-4816-8b88-eaed6987f70f" />



## AFTER (v0.19.0) with Daily Arcs 

<img width="1585" height="1084" alt="image" src="https://github.com/user-attachments/assets/cd03e0af-5c1c-4465-abc9-8097d6b5ee80" />

RMS = √(median residual² + σ²), pooling telluric shifts and velocity-corrected stellar residuals.

| Band      |   N |   Median [Å] |   σ [Å] |   RMS [Å] |
|:----------|----:|-------------:|--------:|----------:|
| u         |  20 |       +0.293 |   0.274 |     0.402 |
| g         |  24 |       +0.223 |   0.278 |     0.357 |
| r         |  24 |       +0.029 |   0.249 |     0.250 |
| i         |  33 |       +0.030 |   0.223 |     0.224 |
| All bands | 101 |       +0.079 |   0.282 |     0.293 |

**Telluric lines in i-band**

<img width="1546" height="446" alt="image" src="https://github.com/user-attachments/assets/2986a5d4-fc8c-402a-be21-3a373085ac71" />

**H-beta stellar line in g-band**

<img width="1541" height="742" alt="image" src="https://github.com/user-attachments/assets/11b422be-dbff-4d55-8ba0-9f3c6aa85802" />


## AFTER (v0.19.0) with Attached Arcs 

<img width="1582" height="1080" alt="image" src="https://github.com/user-attachments/assets/e7bbc269-4b75-47cd-9615-dd58ce7f09e2" />

RMS = √(median residual² + σ²), pooling telluric shifts and velocity-corrected stellar residuals.

| Band      |   N |   Median [Å] |   σ [Å] |   RMS [Å] |
|:----------|----:|-------------:|--------:|----------:|
| u         |  19 |       +0.231 |   0.350 |     0.419 |
| g         |  23 |       +0.144 |   0.395 |     0.421 |
| r         |  23 |       +0.014 |   0.291 |     0.291 |
| i         |  34 |       +0.014 |   0.913 |     0.913 |
| All bands |  99 |       +0.069 |   0.606 |     0.610 |

**Telluric lines in i-band**

<img width="1546" height="442" alt="image" src="https://github.com/user-attachments/assets/9767ef9b-e2c5-4e17-aa43-31bfa319e2ea" />

**H-beta stellar line in g-band**

<img width="1542" height="744" alt="image" src="https://github.com/user-attachments/assets/2e3e194d-0350-44e3-a8cb-56e6359cba22" />

It should be noted here that the wavelength solution, when using the attached arcs, is of sufficient accuracy that attempting to further improve the solution with skylines actually introduces noise and degrades the original accuracy. Therefore, the pipeline now includes a gate so that whenever it recognises an attached arc, it will not attempt to further tune the solution with the skylines. 

## AFTER (v0.19.0) with attached arcs & gate to prevent skyline line adjustments on science frames with attached arcs

<img width="1584" height="1082" alt="image" src="https://github.com/user-attachments/assets/6f572308-fdad-4bb2-8982-3f6b387c16e5" />

RMS = √(median residual² + σ²), pooling telluric shifts and velocity-corrected stellar residuals.

| Band      |   N |   Median [Å] |   σ [Å] |   RMS [Å] |
|:----------|----:|-------------:|--------:|----------:|
| u         |  17 |       +0.279 |   0.252 |     0.376 |
| g         |  25 |       +0.154 |   0.191 |     0.245 |
| r         |  24 |       +0.008 |   0.209 |     0.209 |
| i         |  31 |       -0.061 |   0.191 |     0.201 |
| All bands |  97 |       +0.047 |   0.236 |     0.241 |

**Telluric lines in i-band**

<img width="1546" height="447" alt="image" src="https://github.com/user-attachments/assets/9b6fe572-b98c-4af9-8f08-b0b59f6ea013" />


**H-beta stellar line in g-band**


<img width="1539" height="744" alt="image" src="https://github.com/user-attachments/assets/54eec6f5-dc02-4ed8-a0e4-13eff776ad39" />


