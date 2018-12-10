*This folder has replication files for the paper "Estimating Treatment Effects with Causal Forests: An Application", by Athey and Wager. The paper was written following a workshop organized at the 2018 Atlantic Causal Inference Conference, "Empirical Investigation of Methods for Heterogeneity". Below is a note from the workshop organizers describing the dataset (the note has been lightly edited to remove irrelevant logistical details).*

The attached dataset is based on preliminary data extracted from the National Study of Learning Mindsets (http://mindsetscholarsnetwork.org/about-the-network/current-initatives/national-mindset-study/). This is a study that evalutes a "nudge-like"  intervention to change student behavior at very low cost, but seeks to understand heterogeneity in the intervention. The main goal of the study is to test for heterogeneity in the effect of an intervention designed to instill students with a "growth mindset". A growth mindset is the belief that intelligence can be developed. It is contrasted with a fixed mindset: the belief that intelligence is a fixed trait that is set in stone at birth. (http://mindsetscholarsnetwork.org/learning-mindsets/growth-mindset/). Because growth mindset interventions only teach students to see their school differently, but don't try to change schools themselves, then growth mindset interventions are thought to depend in important ways on the school context. But research has not yet interrogated treatment effect heterogeneity sufficiently. 

The National Study of Learning Mindsets was a randomized controlled trial in a national probability sample of U.S. public high schools. For the purpose of this workshop we created a dataset that instead emulates an observational study with similar characteristics as the National Study (including covariate distributions, data structures, and effect sizes).  Attached is a dataset of about 10,000 students in 76 schools with a simulated outcome Y (a continuous measure of achievement), a binary treatment variable Z indicating receipt of the intervention, and the following 10 covariates at both the student and school level:

- S3 - Students' self-reported expectations for success in the future, a proxy for prior achievement, measured prior to random assignment
- C1 - Categorical variable for student race/ethnicity
- C2 - Categorical variable for student identified gender
- C3 - Categorical variable for student first-generation status (i.e. first in family to go to college)
- XC - School-level categorical variable for urbanicity of the school (i.e. rural, suburban, etc.)
- X1 - School-level mean of students' fixed mindsets, reported prior to random assignment
- X2 - School achievement level, as measured by test scores and college preparation for the previous 4 cohorts of students
- X3 - School racial/ethnic minority composition -- i.e. % black, latino, or native/american
- X4 - School poverty concentration -- i.e. % of students who are from families whose incomes fall below the federal poverty line
- X5 - School size - Total # of students in all four grade levels in the school

The main research questions to be addressed include:

1. Was the mindset intervention effective in improving student achievement?
2. Researchers hypothesize that the effect of the intervention is moderated by school level achievement (X2) and pre-existing mindset norms (X1). In particular there are two competing hypotheses about how X2 moderates the effect of the intervention: Either it is largest in middle-achieving schools (a "Goldilocks effect") or is decreasing in school-level achievement
3. Researchers also collected other covariates and are interested in exploring their possible role in moderating treatment effects.

In your analysis we would ask you to address the the three research questions above. Please recall from our initial invitation that this is not intended to be a "bake off", but rather an opportunity to understand the strengths and weaknesses of methods for addressing important scientific questions, such as those arising in the Mindset Study.  As such we are particularly interested in comments on how you think your chosen methods are either well-suited or require further development to address each of the questions above, rather than specific estimates and/or tests.

Finally, we have worked with the editor of Observational Studies to create a special issue summarizing the findings of this workshop. To that end, all panelists are invited to submit a writeup of their results and comments.

We are looking forward to seeing you all there! This should be fun.

Best,

Carlos Carvalho, Jared Murray, Jennifer Hill and Avi Feller.

NOTE: The dataset may only be used to evaluate methods for inferring treatment effects. Analysis of or publication using other features of the joint distribution of the covariates (the S, C and X variables) are expressly prohibited (and ill adivsed, as noise has been added to the other variables). Any publication using the dataset for its intended purpose should cite the National Study.