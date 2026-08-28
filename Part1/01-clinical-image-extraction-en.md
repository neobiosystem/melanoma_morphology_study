# Part 1 — Clinical Image Feature Extraction (Data Source)

## 1. Positioning and Role in the Closed Loop

This part is the data-source stage of the four-step closed-loop workflow. It supplies the morphological baseline that later stages use for simulation and calibration. Its input is a public dermoscopic and clinical skin-image dataset, and its output is a set of quantitative morphological descriptors, among which the boundary fractal dimension and the compactness index play the central role. These two quantities serve as the calibration targets for the simulation-to-image alignment in Part 3. The scope of this part is confined to objectively extracting descriptors from images; it does not interpret the biological mechanisms behind the descriptors, which is the joint task of Parts 2 and 3. The technical emphasis is therefore on the reproducibility of segmentation and measurement rather than on classification performance as such.

## 2. Method and Technical Route

Image processing is carried out in the MATLAB Image Processing Toolbox, following the sequence of colour-space transformation, lesion segmentation, boundary extraction, and morphological measurement. The raw RGB image is converted to the Lab and HSV spaces, and the luminance or hue channel is used to enhance the contrast between the lesion and the surrounding skin in preparation for thresholding. Segmentation uses an Otsu global threshold to produce an initial mask, with morphological opening and closing to remove artefacts such as hairs and bubbles. For dataset subsets that provide expert-annotated masks, the annotations are used directly to keep the segmentation consistent and reproducible. The lesion boundary is obtained by edge extraction from the segmented mask, using a Canny operator or a morphological gradient, yielding a closed contour.

The fractal dimension is computed on the extracted boundary contour with the box-counting method. The boundary curve is covered with square grids of side length ε, the number N(ε) of boxes intersecting the boundary is counted, and a least-squares line is fitted in the log N(ε) versus log(1/ε) coordinates; the resulting slope is the fractal dimension D_f. The box size ε must be chosen within a scaling region whose lower bound avoids pixel-level noise and whose upper bound avoids the lesion scale, since either truncation would bias the fitted slope. In addition to the fractal dimension, the compactness index (the ratio of perimeter squared to 4π times area), an asymmetry measure, and circularity are computed, forming a set of complementary morphological descriptors so that no single descriptor is over-exposed to one error source.

| Parameter / step | Meaning | Recommended setting |
|---|---|---|
| Colour space | Transformation for contrast enhancement | Lab luminance or HSV hue, depending on image modality |
| Segmentation | Lesion–background separation | Otsu threshold as primary, expert masks as reference, morphological opening/closing for artefacts |
| Boundary extraction | Obtaining a closed contour | Canny operator or morphological gradient |
| Box size ε | Grid size for box-counting | Logarithmically spaced sequence within the scaling region, avoiding pixel-level and lesion-level truncation |
| Fitting | Slope from log N(ε)–log(1/ε) | Least-squares line, recording the coefficient of determination |

## 3. Key Algorithm and Pseudocode

Pseudocode for the box-counting estimate of the boundary fractal dimension follows.

```
Input:  boundary mask B (binary), box-size sequence {e1, ..., eK}
Output: fractal dimension Df, coefficient of determination R2

Normalise B to the [0,1]^2 domain
for k = 1 .. K:
    e = e_k
    N(k) = number of e x e grid boxes intersecting the boundary
    x(k) = log(1/e)
    y(k) = log N(k)
Fit a least-squares line to (x, y)
Df = slope of the fitted line
R2 = coefficient of determination
```

## 4. Expected Output and Evaluation Metrics

The output is the distribution of morphological descriptors for benign nevi and melanomas. Consistent with prior work, benign nevi have smooth boundaries with a fractal dimension close to 1, whereas melanomas show a higher fractal dimension owing to infiltrative growth, so the two classes are expected to be separable in fractal dimension and compactness. This part reports descriptive statistics for the two classes and assesses the discriminative power of each descriptor.

| Metric | Definition | Benign expectation | Malignant expectation |
|---|---|---|---|
| Fractal dimension D_f | Box-counting slope | Close to 1, smooth boundary | Markedly higher, fragmented boundary |
| Compactness | perimeter² / (4π × area) | Close to 1, near-circular | Clearly larger, convoluted contour |
| Asymmetry | Deviation from principal-axis symmetry | Small | Larger |
| Circularity | Shape deviation from an equal-area circle | Higher | Lower |

## 5. Feasibility, Risks, and Alternatives

The fractal dimension and the box-counting method have established precedents in melanoma border analysis, and MATLAB provides complete image-processing support, so the implementation risk is low. The main risk is the effect of segmentation quality on boundary morphology: hairs, bubbles, and uneven illumination in dermoscopic images can contaminate the boundary and perturb the fractal dimension. The alternative is to introduce a hair-removal preprocessing step, to adopt segmentation methods that are more robust to noise such as active contours or graph cuts, and to recompute the descriptors under several segmentation settings for the same lesion to assess robustness. Another risk is the difference in imaging scale between dermoscopic images and clinical photographs, which must be addressed by stating the image modality at the data-selection stage and by normalising the pixel scale before cross-sample comparison.

## 6. Data and Ethical Compliance

This part uses public, de-identified datasets such as ISIC and HAM10000. The images do not contain directly identifiable personal information, so this is secondary use of public data and does not require additional ethics approval, provided the usage terms and citation requirements of each dataset are observed. If the study is later extended to clinical photographs collected at a local institution, ethics-committee review and patient informed consent will be required; the current stage does not involve that scenario.

## 7. References

[1] E. Claridge, P. N. Hall, M. Keefe, and J. P. Allen, "Shape analysis for classification of malignant melanoma," J. Biomed. Eng., vol. 14, no. 3, pp. 229–234, 1992, doi: 10.1016/0141-5425(92)90057-R.

[2] S. H. Carbonetto and S. E. Lew, "Characterization of border structure using fractal dimension in melanomas," in Proc. Annu. Int. Conf. IEEE Eng. Med. Biol. Soc. (EMBC), 2010, doi: 10.1109/IEMBS.2010.5627296.

[3] A.-R. Ali, J. Li, G. Yang, and S. J. O'Shea, "A machine learning approach to automatic detection of irregularity in skin lesion border using dermoscopic images," PeerJ Comput. Sci., vol. 6, e268, 2020, doi: 10.7717/peerj-cs.268.

[4] P. Tschandl, C. Rosendahl, and H. Kittler, "The HAM10000 dataset, a large collection of multi-source dermatoscopic images of common pigmented skin lesions," Sci. Data, vol. 5, 180161, 2018, doi: 10.1038/sdata.2018.161.

[5] N. C. F. Codella et al., "Skin lesion analysis toward melanoma detection 2018: A challenge hosted by the International Skin Imaging Collaboration (ISIC)," arXiv:1902.03368, 2019.
