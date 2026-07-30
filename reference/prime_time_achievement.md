# Indiana Prime Time Third Grade Achievement Evaluation Data

The complete student level data file from the 2000 to 2001 Indiana
Department of Education program evaluation of Project Prime Time,
reported in Lapsley, Daytner, Kelley, and Maxwell (2002, ERIC ED466679).
The evaluation examined the academic performance of *N* = 10,927 third
grade students in *n* = 573 classrooms (here 586 by the
`paste(corp, school, class)` rule), *n* = 163 schools, *n* = 61 school
corporations, and 9 Indiana educational service regions as a function of
class size, pupil to teacher ratio, and the presence of a Prime Time
instructional assistant. The data have been used as a multilevel example
through chapters 3, 4, 6, 9, and 10 of Finch, Bolin, and Kelley (2019,
*Multilevel Modeling Using R*, 2nd ed., CRC Press) and are made
available here as a benchmark data set for the design, measurement, and
analysis of nested data.

## Usage

``` r
prime_time_achievement
```

## Format

A data frame with 10,927 observations on 113 variables. Variables fall
into seven blocks: three derived unique cluster identifiers, student
level demographics and ability and achievement scores, classroom level
variables, school level variables, and school corporation (district)
level variables. Original Indiana DOE variable spellings are preserved,
including the typos `calender` (calendar), `hispanc1` and `hispanc2`
(hispanic), and `rmediate` (remediate), for code compatibility with
Finch, Bolin, and Kelley (2019). The SPSS variable label from the source
file is available as `attr(prime_time_achievement$VAR, "label")` for
every variable carried over from the SPSS file; the one derived recode
without a source label is `classize` (see its entry below).

- `id`:

  Student identifier in the source file (not guaranteed unique; see
  `class_id` etc. for stable cluster keys).

- `region`:

  Indiana educational service region, coded 1 to 9. Sampling stratifier
  (25% of corporations per region).

- `corp`:

  School corporation (district) numeric identifier. **Note:** corp 2400
  appears in both region 2 (12 schools) and region 3 (1 school) in the
  source file, so `paste(region, corp)` is the cleaner cluster key; see
  `corp_id`.

- `school`:

  Numeric school identifier (unique within corporation).

- `class`:

  Classroom number within school, 1 to 8. **Not** unique across schools.

- `corp_id`:

  Derived. `paste(region, corp, sep = "_")`. 61 distinct values,
  matching the count of school corporations reported in Lapsley et al.
  (2002).

- `school_id`:

  Derived. `paste(corp, school, sep = "_")`. 163 distinct values.

- `class_id`:

  Derived. `paste(corp, school, class, sep = "_")`. 586 distinct values.
  (The published report counted 573 classrooms; the small discrepancy
  reflects a different counting convention used in the manuscript.)

- `gender`:

  1 = Female, 2 = Male. 45 `NA`.

- `age`:

  Student age in months.

- `race`:

  Indiana DOE 6-category ethno-racial code: 1 = American Indian /
  Alaskan, 2 = African American, 3 = Asian American, 4 = Hispanic
  American, 5 = Caucasian American, 6 = Multi-racial.

- `geread`:

  Gates-MacGinitie reading.

- `gevocab`:

  Gates-MacGinitie vocabulary.

- `gereadcm`:

  Gates-MacGinitie reading composite.

- `gelang`:

  Gates-MacGinitie language.

- `gelangmc`:

  Gates-MacGinitie language mechanics.

- `gelangcm`:

  Gates-MacGinitie language composite.

- `gemath`:

  Gates-MacGinitie mathematics.

- `gemathcp`:

  Gates-MacGinitie mathematics computation.

- `gemathcm`:

  Gates-MacGinitie mathematics composite.

- `getotal`:

  Gates-MacGinitie total.

- `ncread`:

  NCE reading (ISTEP+).

- `ncvocab`:

  NCE vocabulary (ISTEP+).

- `ncreadcm`:

  NCE reading composite (ISTEP+).

- `nclang`:

  NCE language (ISTEP+).

- `nclangmc`:

  NCE language mechanics (ISTEP+).

- `nclangcm`:

  NCE language composite (ISTEP+).

- `ncmath`:

  NCE mathematics (ISTEP+).

- `ncmathcp`:

  NCE mathematics computation (ISTEP+).

- `ncmathcm`:

  NCE mathematics composite (ISTEP+).

- `nctotal`:

  **NCE total composite (ISTEP+)**, the criterion variable in the
  Lapsley et al. (2002) HLM analyses.

- `aaread`:

  AANCE reading.

- `aavocab`:

  AANCE vocabulary.

- `aareadcm`:

  AANCE reading composite.

- `aalang`:

  AANCE language.

- `aalangmc`:

  AANCE language mechanics.

- `aalangcm`:

  AANCE language composite.

- `aamath`:

  AANCE mathematics.

- `aamathcp`:

  AANCE mathematics computation.

- `aamathcm`:

  AANCE mathematics composite.

- `aatotal`:

  AANCE total.

- `npanverb`:

  NPA nonverbal reasoning.

- `npamem`:

  NPA working memory.

- `npaverb`:

  NPA verbal reasoning.

- `npatotal`:

  NPA total.

- `csi`:

  Cognitive Skills Index (student level).

- `multi`:

  Multi-age classroom indicator (1 = Yes, 2 = No).

- `typmulti`:

  Type of multi-age classroom (1 = 1st-2nd-3rd grades, 2 = 2nd-3rd, 3 =
  3rd-4th, 4 = 2nd-3rd- 4th; `NA` when `multi == 2`).

- `clenroll`:

  Official class enrollment.

- `classize`:

  Project STAR style class size category: 1 = small (roughly 12-17), 2 =
  regular (roughly 18-22), 3 = regular-larger (roughly 23-26), 4 = large
  (27 or more). Boundaries follow the STAR classification (Pate-Bain and
  Achilles, 1986).

- `ptratio`:

  Classroom pupil to teacher ratio (IDOE formula: enrollment / \[1.00
  per full time teacher + 0.33 per full time aide + 0.165 per part time
  aide\]).

- `ptia`:

  Prime Time Instructional Aide status: 1 = aide present, 2 = no aide, 3
  = other assistant listed. The focal treatment indicator.

- `ptstatus`:

  Status of Prime Time aide: 1 = full time in classroom, 2 = part time
  in classroom; `NA` when `ptia != 1`.

- `locale`:

  NCES locale code (1 = large central city, 2 = mid-size central city, 3
  = urban fringe of large city, 4 = urban fringe of mid-size city, 5 =
  large town, 6 = small town, 7 = rural).

- `chapter1`:

  School receives Title I (legacy "Chapter 1") money? 1 = Yes, 2 = No.

- `ses`:

  School SES, the IDOE percentage of students *not* eligible for
  subsidized lunch (0 to 100; higher = more affluent). The school level
  SES variable used in the Lapsley et al. (2002) HLM analyses.

- `context`:

  IDOE contextual rank for the school.

- `calender`:

  School calendar type (1 = traditional, 2 = year round). Typo preserved
  from source.

- `senroll`:

  Building (school) enrollment.

- `sattend`:

  Building attendance rate (percent).

- `white1`:

  School percent White.

- `black1`:

  School percent Black.

- `hispanc1`:

  School percent Hispanic (typo preserved).

- `asian1`:

  School percent Asian.

- `aindian1`:

  School percent American Indian.

- `multi1`:

  School percent multi-racial.

- `total1`:

  School total percent non-white.

- `noteach`:

  Number of teachers in the building (full time equivalent).

- `avgage1`:

  School average teacher age.

- `avgexp1`:

  School average teacher experience (years).

- `avgsal1`:

  School average teacher salary (dollars).

- `spert`:

  School students per teacher.

- `thrdclss`:

  Number of third grade classrooms in the building.

- `thrdstud`:

  Number of third graders who took ISTEP+ in the building.

- `passla1`:

  Building percent passing language arts.

- `passmth1`:

  Building percent passing math.

- `passbth1`:

  Building percent passing both.

- `tmnnce1`:

  Building total battery mean NCE.

- `rmdnce1`:

  Building reading *median* NCE. **Note:** the 56 rows from one building
  (`corp` 5740, `school` 6187) carry the source value 7603, an evident
  data-entry error in the Indiana DOE file (an NCE is on the 1 to 99
  scale, and the same building's other median-NCE columns are in range).
  The value is preserved as shipped rather than silently corrected,
  since the true value cannot be recovered; drop or set it to `NA`
  before analyzing this column.

- `lamdnce1`:

  Building language arts *median* NCE.

- `mmdnce1`:

  Building mathematics *median* NCE.

- `tmdnce1`:

  Building total battery *median* NCE.

- `avgcsi1`:

  Building average Cognitive Skills Index.

- `geog`:

  Geographic category of the corporation (1 = urban, 2 = suburban, 3 =
  town, 4 = rural). Sampling stratifier within region.

- `totepp`:

  Corporation total expense per pupil (1997–1999, dollars).

- `cenroll`:

  Corporation enrollment (all grades).

- `cattend`:

  Corporation attendance rate (percent).

- `freelnch`:

  Corporation percent eligible for free lunch.

- `lep`:

  Corporation percent with limited English proficiency.

- `speced`:

  Corporation percent in special education.

- `minority`:

  Corporation percent minority.

- `white2`:

  Corporation total White public enrollment (raw count).

- `black2`:

  Corporation total Black public enrollment.

- `hispanc2`:

  Corporation total Hispanic public enrollment (typo preserved).

- `asian2`:

  Corporation total Asian public enrollment.

- `aindian2`:

  Corporation total American Indian public enrollment.

- `multi2`:

  Corporation total multi-racial public enrollment.

- `total2`:

  Corporation total non-white public enrollment.

- `thrdadm`:

  Corporation third grade ADM (average daily membership).

- `thrdtech`:

  Corporation third grade teachers.

- `avgage2`:

  Corporation average teacher age.

- `avgexp2`:

  Corporation average teacher experience (years).

- `avgsal2`:

  Corporation average teacher salary (dollars).

- `thrdaide`:

  Corporation third grade aides.

- `passla2`:

  Corporation percent passing language arts.

- `passmth2`:

  Corporation percent passing math.

- `passbth2`:

  Corporation percent passing both.

- `tmnnce2`:

  Corporation total battery mean NCE.

- `rmdnce2`:

  Corporation reading *median* NCE.

- `lamdnce2`:

  Corporation language arts *median* NCE.

- `mmdnce2`:

  Corporation mathematics *median* NCE.

- `tmdnce2`:

  Corporation total battery *median* NCE.

- `rmediate`:

  Corporation remediation funding per pupil (dollars). Typo preserved.

## Source

Indiana Department of Education program evaluation of Project Prime
Time, 2000 to 2001 academic year. Sample of 10,927 third grade students
in 586 classrooms in 163 schools in 61 corporations in 9 educational
service regions. The records are public data that the author, a member
of the evaluation team, is authorized to distribute.

Lapsley, D. K., Daytner, K. M., Kelley, K., and Maxwell, S. E. (2002).
*Teacher aides, class size and academic achievement: A preliminary
evaluation of Indiana's Prime Time*. Paper presented at the Annual
Meeting of the American Educational Research Association, New Orleans,
LA, April 1-5, 2002. ERIC document ED466679.

## Details

The data frame is exposed under two names for user convenience:
`prime_time_achievement` (canonical) and the short alias `Prime_Time`.
The canonical name is the one shipped in `data/`, so
`data(prime_time_achievement)` loads it; the alias `Prime_Time` is an
active binding created at load time (rather than a duplicate `.rda`), so
the bare name `Prime_Time` resolves to the identical object after
[`library(DMAR)`](https://kenkelley.org) but `data(Prime_Time)` does not
find a data set. A parallel lowercase `prime_time` alias is
intentionally *not* provided because shipping both `prime_time.rda` and
`Prime_Time.rda` would collide on case-insensitive filesystems (APFS,
HFS+, NTFS).

**Study background.** Indiana's Prime Time program, phased in beginning
1984 to 1985 (Indiana statute; House Bill 1166 of 2001 codified the
modern funding formula), was one of the earliest state level initiatives
in the United States to use a funding formula to reduce class size and
pupil to teacher ratio in Kindergarten through third grade. Funds were
distributed to school corporations to maintain a corporation average
pupil to teacher ratio of 18:1 in K and grade 1 and 20:1 in grades 2 and
3; corporations could meet the target by hiring additional teachers or,
more commonly, paraprofessional instructional assistants. Along with
Tennessee's Project STAR (Pate-Bain and Achilles, 1986; covered by
Education Week, *Research: Sizing Up Small Classes*, February 2001),
Prime Time was a widely cited national model.

In 1999 the Indiana Department of Education funded a three year program
evaluation of Prime Time. The third year of the evaluation, conducted by
Daniel K. Lapsley and Katrina M. Daytner (Ball State University and
Western Illinois University) with technical assistance from Ken Kelley
and Scott E. Maxwell (University of Notre Dame), examined the academic
performance of randomly selected Indiana third graders on the state
mandated ISTEP+ standardized achievement test as a function of class
size, pupil to teacher ratio, and the presence of a Prime Time
instructional aide, using hierarchical linear modeling. The preliminary
report and the AERA 2002 paper that summarizes the analyses are archived
as ERIC document ED466679 (Lapsley, Daytner, Kelley, and Maxwell, 2002).
An earlier background report from the same evaluation team, prior to the
Notre Dame group joining the project, is archived as ERIC document
ED455220.

**Sampling.** School corporations were drawn by stratified cluster
sampling with two rules: 25% of corporations from each of the nine
Indiana educational service regions, and at least one urban corporation
per region, with the remainder proportionally allocated across
geographic categories (urban, suburban, town, rural). The achieved
sample was 61 corporations (78% of the target), 163 schools, 573
classrooms as counted in the manuscript (586 by the
`paste(corp, school, class)` rule used here), and 10,927 students (49.6%
female; 85% Caucasian, 9.2% African American, 3.2% Hispanic). 4,016
students were in classrooms with a Prime Time instructional assistant
(`ptia == 1`), 6,765 in classrooms without (`ptia == 2`); the file here
shows 4,021 and 6,789 plus 117 `ptia == 3` (other assistant listed),
with the small differences reflecting cleaning rules applied between the
manuscript count and the final SPSS file.

**Instruments.** Third graders sit for the ISTEP+ (Indiana Statewide
Testing for Educational Progress) in September of the school term. The
ISTEP+ is published by CTB/McGraw-Hill and includes language arts,
reading, and mathematics assessments. Normal Curve Equivalent (NCE)
composite scores for these domains and for the total are the
`ncread / nclang / ncmath / nctotal` columns and were the criterion
variables in the published HLM analyses. NCE scores have a population
mean of 50 and a standard deviation of approximately 21.06, with
percentiles 1, 50, and 99 mapping to NCE scores of 1, 50, and 99. The
`ge*` family is the parallel Gates-MacGinitie battery; the `aa*` family
is the African American comparison NCE (AANCE); the `npa*` family is the
cognitive abilities battery used as student level covariates in Finch,
Bolin, and Kelley (2019).

**Nested data structure.** The natural hierarchy is student *within*
classroom *within* school *within* corporation *within* region. The
derived identifiers `corp_id`, `school_id`, and `class_id` are
pre-computed and safe to use as grouping variables; the bare `corp` and
`class` columns are *not* unique by themselves. Class sizes range from 3
to 28 students (median 19); schools have 1 to 8 third grade classrooms
(median 3) and 11 to 166 students (median 65); corporations have 15 to
808 students (median 124). Variance decomposition for the published
outcome `nctotal` based on the three level random intercept null model
`lmer(nctotal ~ 1 + (1 | corp_id/school_id))` gives: between-corporation
variance 16.50, between-school within corporation variance 22.70, and
within school residual variance 240.43, so that
*ICC*\\\_{\mathrm{corp}}\\ \\\approx\\ 0.059,
*ICC*\\\_{\mathrm{school\|corp}}\\ \\\approx\\ 0.081, and the combined
cluster *ICC*\\\_{\mathrm{cluster}}\\ \\\approx\\ 0.140. These
nontrivial intraclass correlations are the methodological reason
multilevel modeling is preferred to ordinary least squares regression
for these data.

**Level 1, 2, 3 model framework.** For an outcome \\Y\_{ijk}\\ on
student \\i\\ in classroom \\j\\ in school \\k\\, with student level
predictor \\X^{(1)}\_{ijk}\\, classroom level predictor
\\X^{(2)}\_{jk}\\, and school level predictor \\X^{(3)}\_k\\, the
published Lapsley et al. (2002) family of HLM models has the equations
\$\$Y\_{ijk} = \pi\_{0jk} + \pi\_{1jk} X^{(1)}\_{ijk} + e\_{ijk} \quad
\text{(Level 1)},\$\$ \$\$\pi\_{0jk} = \beta\_{00k} + \beta\_{01k}
X^{(2)}\_{jk} + r\_{0jk}, \quad \pi\_{1jk} = \beta\_{10k} + r\_{1jk}
\quad \text{(Level 2)},\$\$ \$\$\beta\_{00k} = \gamma\_{000} +
\gamma\_{001} X^{(3)}\_k + u\_{00k}, \quad \beta\_{01k} = \gamma\_{010},
\quad \beta\_{10k} = \gamma\_{100} \quad \text{(Level 3)},\$\$ with
\\e\_{ijk} \sim N(0, \sigma^2)\\, \\r\_{jk} \sim N(0,
\mathbf{T}\_\pi)\\, and \\u\_{00k} \sim N(0, \tau\_{00})\\. Substituting
upward, the reduced form is \$\$Y\_{ijk} = \gamma\_{000} + \gamma\_{100}
X^{(1)}\_{ijk} + \gamma\_{010} X^{(2)}\_{jk} + \gamma\_{001}
X^{(3)}\_k + u\_{00k} + r\_{0jk} + r\_{1jk} X^{(1)}\_{ijk} +
e\_{ijk},\$\$ which in lme4 translates to
`lmer(Y ~ X1 + X2 + X3 + (1 + X1 | corp_id/school_id))`. See the
examples for concrete fits.

**Suggested benchmark uses.** The data set is intentionally rich enough
to support a wide range of demonstrations and benchmarks, including:

- Two-, three-, and four-level random intercept and random slope models
  with lme4, nlme, or glmmTMB.

- Cross-level interaction modeling (e.g., race \\\times\\ class size,
  ptia \\\times\\ SES).

- ICC, design effect, and cluster level sample size calculations.

- Comparisons of unweighted vs. design weighted estimators for
  stratified cluster samples.

- Bayesian multilevel modeling and prior sensitivity (brms, MCMCglmm,
  rstanarm).

- Missing data demonstrations (the `ge*`, `nc*`, and `aa*` columns have
  non-trivial missingness; see
  `vapply(prime_time_achievement, function(x) sum(is.na(x)), integer(1))`).

- Multilevel reliability, intraclass correlation, and measurement
  invariance demonstrations across schools, corporations, and the
  categorical predictors.

**Privacy and identifiability.** The student level rows contain no
names, addresses, or other personally identifiable information.
Demographic variables are age in months, gender, and a six category race
code; all other fields are test scores or aggregated school /
corporation statistics. The numeric `corp`, `school`, and `class`
identifiers are the same administrative numbers used in the original
Indiana Department of Education public files for the 2000 to 2001 school
year; they could in principle be cross referenced to that public
information to identify specific schools or corporations. No individual
student can be identified from any combination of variables in this
file.

**Missing data convention.** The Indiana DOE source used `999` as the
student level missing data code and `888` as the "not applicable" code
for `typmulti` and `ptstatus`. The build script converts both to `NA`
(the SPSS missingness ranges already do most of the recoding on import).
The retained SPSS variable label is available via
`attr(prime_time_achievement$X, "label")` on every variable carried over
from the SPSS file (all columns except the derived recode `classize`).

## References

*Primary citation.* Lapsley, D. K., Daytner, K. M., Kelley, K., and
Maxwell, S. E. (2002). *Teacher aides, class size and academic
achievement: A preliminary evaluation of Indiana's Prime Time*. ERIC
document ED466679. <https://eric.ed.gov/?id=ED466679>.

*Use as a multilevel modeling running example.* Finch, W. H., Bolin, J.
E., and Kelley, K. (2019). *Multilevel modeling using R* (2nd ed.). CRC
Press. The 2nd edition (Finch, Bolin, and Kelley, 2019) is the edition
that uses these data; later editions are not authored by Kelley and
should not be cited for that use.

*Background report from the same evaluation team.* Lapsley, D. K., and
Daytner, K. M. (2001). Indiana's class size reduction initiative:
Teacher perspectives on training, implementation, and pedagogy. ERIC
document ED455220. <https://files.eric.ed.gov/fulltext/ED455220.pdf>.

*Indiana statutory context.* Indiana General Assembly, House Bill 1166
(2001). <https://archive.iga.in.gov/2001/bills/IN/IN1166.1.html>.

*Project STAR background and Education Week coverage.* Pate-Bain, H.,
and Achilles, C. M. (1986). Interesting developments on class size. *Phi
Delta Kappan, 67*, 662–665. See also Education Week, *Research: Sizing
up small classes* (February 7, 2001),
<https://www.edweek.org/leadership/research-sizing-up-small-classes/2001/02>.

*Project STAR teacher aide null result that motivated the Prime Time
evaluation.* Finn, J. D., Gerber, S. B., Farber, S. L., and Achilles, C.
M. (2000). Teacher aides: An alternative to small classes? In M. C. Wang
and J. D. Finn (Eds.), *How small classes help teachers do their best*
(pp. 131–174). Temple University Center for Research in Human
Development and Education.

## Author

Ken Kelley

## Examples

``` r
data(prime_time_achievement)
dim(prime_time_achievement)
#> [1] 10927   113

# Variable labels from the SPSS source are preserved on every column:
attr(prime_time_achievement$nctotal, "label")
#> [1] "NCE TOTAL"
attr(prime_time_achievement$ptia,    "label")
#> [1] "PRESENCE OF A PRIME TIME IA?"

# Cluster counts (reconciled with Lapsley et al., 2002):
length(unique(prime_time_achievement$corp_id))    # 61
#> [1] 61
length(unique(prime_time_achievement$school_id))  # 163
#> [1] 163
length(unique(prime_time_achievement$class_id))   # 586
#> [1] 586

# Reconciling with the manuscript:
table(prime_time_achievement$gender, useNA = "ifany")
#> 
#>    1    2 <NA> 
#> 5425 5457   45 
table(prime_time_achievement$race,   useNA = "ifany")
#> 
#>    1    2    3    4    5    6 <NA> 
#>   16  995   62  348 9207  188  111 
table(prime_time_achievement$ptia)
#> 
#>    1    2    3 
#> 4021 6789  117 
table(prime_time_achievement$classize)
#> 
#>    1    2    3    4 
#> 1085 4571 4854  417 

# ----- Selecting subsets of interest -----

# Caucasian and African American only (the matched-race
# supplementary analyses in Lapsley et al., 2002):
pt_wb <- subset(prime_time_achievement, race %in% c(2L, 5L))

# Drop the few "other assistant listed" cases for a clean
# aide / no-aide contrast:
pt_clean <- subset(prime_time_achievement, ptia %in% c(1L, 2L))

# Only rural corporations (geog == 4), which is what the source
# SPSS file's FILTER_$ variable encoded:
pt_rural <- subset(prime_time_achievement, geog == 4L)

# Complete cases on the nctotal-on-race-and-class-size analysis:
analysis_vars <- c("nctotal", "race", "classize", "ses",
                   "corp_id", "school_id", "class_id")
pt_complete <- prime_time_achievement[
  complete.cases(prime_time_achievement[, analysis_vars]),
  analysis_vars
]

# \donttest{
# Three-level null random intercept model. The variance
# decomposition gives the corp, school | corp, and within
# ICCs reported in the Details section.
m_null <- lme4::lmer(nctotal ~ 1 + (1 | corp_id/school_id),
                     data = prime_time_achievement)
summary(m_null)
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: nctotal ~ 1 + (1 | corp_id/school_id)
#>    Data: prime_time_achievement
#> 
#> REML criterion at convergence: 90763.2
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.6512 -0.7367 -0.0133  0.7201  2.8518 
#> 
#> Random effects:
#>  Groups            Name        Variance Std.Dev.
#>  school_id:corp_id (Intercept)  22.72    4.766  
#>  corp_id           (Intercept)  16.29    4.037  
#>  Residual                      240.43   15.506  
#> Number of obs: 10865, groups:  school_id:corp_id, 163; corp_id, 61
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept)  60.5874     0.7061    85.8

# Level-1 (race), level-2 (ptia and classize), and level-3
# (ses) main-effect model. Compare to Lapsley et al. (2002),
# which fit closely related HLM specifications.
m_main <- lme4::lmer(
  nctotal ~ factor(race) + factor(ptia) + classize + ses +
    (1 | corp_id/school_id),
  data = prime_time_achievement)
summary(m_main)
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: nctotal ~ factor(race) + factor(ptia) + classize + ses + (1 |  
#>     corp_id/school_id)
#>    Data: prime_time_achievement
#> 
#> REML criterion at convergence: 89620.2
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.4991 -0.7372 -0.0115  0.7199  2.9917 
#> 
#> Random effects:
#>  Groups            Name        Variance Std.Dev.
#>  school_id:corp_id (Intercept)  16.24    4.029  
#>  corp_id           (Intercept)  13.61    3.689  
#>  Residual                      236.75   15.387  
#> Number of obs: 10755, groups:  school_id:corp_id, 163; corp_id, 61
#> 
#> Fixed effects:
#>               Estimate Std. Error t value
#> (Intercept)   49.24631    4.49485  10.956
#> factor(race)2 -8.26722    4.02948  -2.052
#> factor(race)3  7.36730    4.45236   1.655
#> factor(race)4 -3.67271    4.12163  -0.891
#> factor(race)5 -0.30283    3.98688  -0.076
#> factor(race)6 -1.99673    4.14749  -0.481
#> factor(ptia)2  0.64273    0.66173   0.971
#> factor(ptia)3 -0.56790    2.64478  -0.215
#> classize       1.46167    0.35583   4.108
#> ses            0.11198    0.02294   4.881
#> 
#> Correlation of Fixed Effects:
#>             (Intr) fctr(r)2 fctr(r)3 fct()4 fct()5 fct()6 fctr(p)2 fctr(p)3
#> factor(rc)2 -0.885                                                         
#> factor(rc)3 -0.791  0.883                                                  
#> factor(rc)4 -0.863  0.963    0.864                                         
#> factor(rc)5 -0.885  0.986    0.893    0.965                                
#> factor(rc)6 -0.855  0.953    0.860    0.934  0.959                         
#> factor(pt)2 -0.183  0.000   -0.006   -0.005 -0.002 -0.004                  
#> factor(pt)3 -0.050  0.001   -0.001    0.001 -0.002  0.000  0.126           
#> classize    -0.183  0.009    0.011    0.000  0.006  0.001  0.265    0.011  
#> ses         -0.377  0.023   -0.004    0.018 -0.004  0.010  0.112    0.071  
#>             classz
#> factor(rc)2       
#> factor(rc)3       
#> factor(rc)4       
#> factor(rc)5       
#> factor(rc)6       
#> factor(pt)2       
#> factor(pt)3       
#> classize          
#> ses         -0.090

# Cross-level interaction: ptia x ses (the published finding
# was that aide benefit was concentrated in higher-SES
# schools).
m_inter <- lme4::lmer(
  nctotal ~ factor(race) + factor(ptia) * ses + classize +
    (1 | corp_id/school_id),
  data = prime_time_achievement)
summary(m_inter)
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: nctotal ~ factor(race) + factor(ptia) * ses + classize + (1 |  
#>     corp_id/school_id)
#>    Data: prime_time_achievement
#> 
#> REML criterion at convergence: 89623.4
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.4958 -0.7366 -0.0118  0.7197  2.9910 
#> 
#> Random effects:
#>  Groups            Name        Variance Std.Dev.
#>  school_id:corp_id (Intercept)  16.39    4.049  
#>  corp_id           (Intercept)  13.58    3.685  
#>  Residual                      236.76   15.387  
#> Number of obs: 10755, groups:  school_id:corp_id, 163; corp_id, 61
#> 
#> Fixed effects:
#>                    Estimate Std. Error t value
#> (Intercept)       48.659313   5.368981   9.063
#> factor(race)2     -8.266518   4.029644  -2.051
#> factor(race)3      7.359468   4.452546   1.653
#> factor(race)4     -3.673089   4.121863  -0.891
#> factor(race)5     -0.307672   3.987043  -0.077
#> factor(race)6     -2.000371   4.147672  -0.482
#> factor(ptia)2      1.259289   3.624373   0.347
#> factor(ptia)3     27.518709  37.753208   0.729
#> ses                0.118886   0.042988   2.766
#> classize           1.461474   0.356267   4.102
#> factor(ptia)2:ses -0.007379   0.044300  -0.167
#> factor(ptia)3:ses -0.423564   0.568556  -0.745
#> 
#> Correlation of Fixed Effects:
#>             (Intr) fctr(r)2 fctr(r)3 fct()4 fct()5 fct()6 fctr(p)2 fctr(p)3
#> factor(rc)2 -0.741                                                         
#> factor(rc)3 -0.662  0.883                                                  
#> factor(rc)4 -0.720  0.963    0.864                                         
#> factor(rc)5 -0.740  0.986    0.893    0.965                                
#> factor(rc)6 -0.714  0.953    0.860    0.934  0.959                         
#> factor(pt)2 -0.565  0.000   -0.003   -0.005 -0.002 -0.004                  
#> factor(pt)3 -0.085  0.001   -0.001    0.001 -0.001  0.000  0.110           
#> ses         -0.631  0.012   -0.004    0.006 -0.003  0.003  0.841    0.121  
#> classize    -0.138  0.009    0.011    0.000  0.006  0.001  0.023   -0.010  
#> fctr(pt)2:s  0.546  0.000    0.002    0.005  0.002  0.003 -0.983   -0.103  
#> fctr(pt)3:s  0.077 -0.001    0.001   -0.001  0.001  0.000 -0.099   -0.997  
#>             ses    classz fc()2:
#> factor(rc)2                     
#> factor(rc)3                     
#> factor(rc)4                     
#> factor(rc)5                     
#> factor(rc)6                     
#> factor(pt)2                     
#> factor(pt)3                     
#> ses                             
#> classize    -0.070              
#> fctr(pt)2:s -0.844  0.026       
#> fctr(pt)3:s -0.110  0.010  0.093
# }
```
