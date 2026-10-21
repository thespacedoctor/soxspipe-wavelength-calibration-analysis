# SOXS UV-VIS Wavelength Calibration Analysis

---

**TL;DR:**

* With daily arcs, we achieve a UV-VIS wavelength solution accuracy of ±0.5Å (for the bluest regions, with greater accuracy achieved as wavelength increases)
* With attached arcs, we achieve a UV-VIS wavelength solution accuracy of ±0.4Å (for the bluest regions, with greater accuracy achieved as wavelength increases)

---

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

<img width="1583" height="1080" alt="image" src="https://github.com/user-attachments/assets/1e86e5fc-3c2c-422b-b982-13f2b5c897d0" />

RMS = √(median residual² + σ²), pooling telluric shifts and velocity-corrected stellar residuals.

| Band      |   N |   Median [Å] |   σ [Å] |   RMS [Å] |
|:----------|----:|-------------:|--------:|----------:|
| u         |  14 |       -0.175 |   0.380 |     0.419 |
| g         |  25 |       -0.418 |   0.618 |     0.746 |
| r         |  21 |       +0.184 |   0.398 |     0.438 |
| i         |  31 |       +0.202 |   0.612 |     0.644 |
| All bands |  91 |       +0.046 |   0.634 |     0.636 |

**Telluric lines in i-band**

<img width="1548" height="441" alt="image" src="https://github.com/user-attachments/assets/e441858a-4547-4ebd-9aac-fb7ec627fc94" />

**H-beta stellar line in g-band**

<img width="1542" height="745" alt="image" src="https://github.com/user-attachments/assets/fab0897d-fd65-460b-968a-7d7559aad760" />


## BEFORE (v0.17.4) with Attached Arcs 

<img width="1583" height="1080" alt="image" src="https://github.com/user-attachments/assets/fce45159-ec8c-4fe8-8e5e-1d1cb5c51f88" />


RMS = √(median residual² + σ²), pooling telluric shifts and velocity-corrected stellar residuals.

| Band      |   N |   Median [Å] |   σ [Å] |   RMS [Å] |
|:----------|----:|-------------:|--------:|----------:|
| u         |  16 |       +0.206 |   0.388 |     0.439 |
| g         |  25 |       +0.017 |   0.678 |     0.678 |
| r         |  22 |       +0.021 |   0.206 |     0.207 |
| i         |  29 |       -0.084 |   0.226 |     0.241 |
| All bands |  92 |       +0.016 |   0.433 |     0.434 |

**Telluric lines in i-band**

<img width="1534" height="749" alt="image" src="https://github.com/user-attachments/assets/dfab0436-c115-47a2-8158-4d2e7b6563ae" />

**H-beta stellar line in g-band**

<img width="1543" height="447" alt="image" src="https://github.com/user-attachments/assets/2da77ddd-4ecb-4816-8b88-eaed6987f70f" />



## AFTER (v0.19.0) with Daily Arcs 

<img width="1585" height="1081" alt="image" src="https://github.com/user-attachments/assets/2b4a78db-4b87-4f5b-a674-b31db810d3cd" />

RMS = √(median residual² + σ²), pooling telluric shifts and velocity-corrected stellar residuals.

| Band      |   N |   Median [Å] |   σ [Å] |   RMS [Å] |
|:----------|----:|-------------:|--------:|----------:|
| u         |  19 |       +0.312 |   0.282 |     **0.420** |
| g         |  24 |       +0.223 |   0.278 |     0.357 |
| r         |  24 |       +0.029 |   0.249 |     0.250 |
| i         |  32 |       +0.026 |   0.216 |     0.218 |
| All bands |  99 |       +0.070 |   0.283 |     0.291 |

**Telluric lines in i-band**

<img width="1546" height="446" alt="image" src="https://github.com/user-attachments/assets/2986a5d4-fc8c-402a-be21-3a373085ac71" />

**H-beta stellar line in g-band**

<img width="1541" height="742" alt="image" src="https://github.com/user-attachments/assets/11b422be-dbff-4d55-8ba0-9f3c6aa85802" />


## AFTER (v0.19.0) with Attached Arcs 

<img width="1583" height="1085" alt="image" src="https://github.com/user-attachments/assets/17ec5f57-de4b-4c95-88a4-1f1da4a53709" />

RMS = √(median residual² + σ²), pooling telluric shifts and velocity-corrected stellar residuals.

| Band      |   N |   Median [Å] |   σ [Å] |   RMS [Å] |
|:----------|----:|-------------:|--------:|----------:|
| u         |  17 |       +0.230 |   0.371 |     0.437 |
| g         |  23 |       +0.142 |   0.399 |     0.423 |
| r         |  22 |       +0.025 |   0.294 |     0.295 |
| i         |  32 |       +0.010 |   0.252 |     0.252 |
| All bands |  94 |       +0.053 |   0.334 |     0.339 |

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
| u         |  17 |       +0.279 |   0.252 |     **0.376** |
| g         |  25 |       +0.154 |   0.191 |     0.245 |
| r         |  24 |       +0.008 |   0.209 |     0.209 |
| i         |  31 |       -0.061 |   0.191 |     0.201 |
| All bands |  97 |       +0.047 |   0.236 |     0.241 |

**Telluric lines in i-band**

<img width="1546" height="447" alt="image" src="https://github.com/user-attachments/assets/9b6fe572-b98c-4af9-8f08-b0b59f6ea013" />


**H-beta stellar line in g-band**


<img width="1539" height="744" alt="image" src="https://github.com/user-attachments/assets/54eec6f5-dc02-4ed8-a0e4-13eff776ad39" />


