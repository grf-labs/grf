# This file reads in replication data from
# Carvalho, L. S., Meier, S., & Wang, S. W. (2016). Poverty and economic decision-making: Evidence from changes in financial resources at payday. American economic review, 106(2), 260-84.
# Farbmacher, H., Kögel, H., & Spindler, M. (2021). Heterogeneous effects of poverty on attention. Labour Economics, 71, 102028.
# and stores processed covariate/treatment/outcome info in carvalho2016.csv

# Data license copy/paste:

# Modified BSD License (https://opensource.org/licenses/BSD-3-Clause)
# - applies to all code, scripts, programs, and SOFTWARE.
# This is any statements or instructions to be used directly or
# indirectly in a computer in order to bring about a certain result,
# and may include interpretive, object or source code.
#
# Creative Commons Attribution 4.0 International Public License
# (https://creativecommons.org/licenses/by/4.0/)
# - applies to databases, images, tables, text, and any other objects
#
# COPYRIGHT 2016 American Economic Association
#
# =================================================================
#   Modified BSD License
# =================================================================
#
#   Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
#   1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
# =================================================================
#   Creative Commons Attribution 4.0 International Public License
# =================================================================
#
#   By exercising the Licensed Rights (defined below), You accept and agree to be
# bound by the terms and conditions of this Creative Commons Attribution 4.0
# International Public License ("Public License"). To the extent this Public
# License may be interpreted as a contract, You are granted the Licensed Rights in
# consideration of Your acceptance of these terms and conditions, and the Licensor
# grants You such rights in consideration of benefits the Licensor receives from
# making the Licensed Material available under these terms and conditions.
#
# Section 1 – Definitions.
#
# Adapted Material means material subject to Copyright and Similar Rights that is
# derived from or based upon the Licensed Material and in which the Licensed
# Material is translated, altered, arranged, transformed, or otherwise modified in
# a manner requiring permission under the Copyright and Similar Rights held by the
# Licensor. For purposes of this Public License, where the Licensed Material is a
# musical work, performance, or sound recording, Adapted Material is always
# produced where the Licensed Material is synched in timed relation with a moving
# image. Adapter's License means the license You apply to Your Copyright and
# Similar Rights in Your contributions to Adapted Material in accordance with the
# terms and conditions of this Public License. Copyright and Similar Rights means
# copyright and/or similar rights closely related to copyright including, without
# limitation, performance, broadcast, sound recording, and Sui Generis Database
# Rights, without regard to how the rights are labeled or categorized. For
# purposes of this Public License, the rights specified in Section 2(b)(1)-(2) are
# not Copyright and Similar Rights. Effective Technological Measures means those
# measures that, in the absence of proper authority, may not be circumvented under
# laws fulfilling obligations under Article 11 of the WIPO Copyright Treaty
# adopted on December 20, 1996, and/or similar international agreements.
# Exceptions and Limitations means fair use, fair dealing, and/or any other
# exception or limitation to Copyright and Similar Rights that applies to Your use
# of the Licensed Material. Licensed Material means the artistic or literary work,
# database, or other material to which the Licensor applied this Public License.
# Licensed Rights means the rights granted to You subject to the terms and
# conditions of this Public License, which are limited to all Copyright and
# Similar Rights that apply to Your use of the Licensed Material and that the
# Licensor has authority to license. Licensor means the individual(s) or
# entity(ies) granting rights under this Public License. Share means to provide
# material to the public by any means or process that requires permission under
# the Licensed Rights, such as reproduction, public display, public performance,
# distribution, dissemination, communication, or importation, and to make material
# available to the public including in ways that members of the public may access
# the material from a place and at a time individually chosen by them. Sui Generis
# Database Rights means rights other than copyright resulting from Directive
# 96/9/EC of the European Parliament and of the Council of 11 March 1996 on the
# legal protection of databases, as amended and/or succeeded, as well as other
# essentially equivalent rights anywhere in the world. You means the individual or
# entity exercising the Licensed Rights under this Public License. Your has a
# corresponding meaning. Section 2 – Scope.
#
# License grant. Subject to the terms and conditions of this Public License, the
# Licensor hereby grants You a worldwide, royalty-free, non-sublicensable,
# non-exclusive, irrevocable license to exercise the Licensed Rights in the
# Licensed Material to: reproduce and Share the Licensed Material, in whole or in
# part; and produce, reproduce, and Share Adapted Material. Exceptions and
# Limitations. For the avoidance of doubt, where Exceptions and Limitations apply
# to Your use, this Public License does not apply, and You do not need to comply
# with its terms and conditions. Term. The term of this Public License is
# specified in Section 6(a). Media and formats; technical modifications allowed.
# The Licensor authorizes You to exercise the Licensed Rights in all media and
# formats whether now known or hereafter created, and to make technical
# modifications necessary to do so. The Licensor waives and/or agrees not to
# assert any right or authority to forbid You from making technical modifications
# necessary to exercise the Licensed Rights, including technical modifications
# necessary to circumvent Effective Technological Measures. For purposes of this
# Public License, simply making modifications authorized by this Section 2(a)(4)
# never produces Adapted Material. Downstream recipients. Offer from the Licensor
# – Licensed Material. Every recipient of the Licensed Material automatically
# receives an offer from the Licensor to exercise the Licensed Rights under the
# terms and conditions of this Public License. No downstream restrictions. You may
# not offer or impose any additional or different terms or conditions on, or apply
# any Effective Technological Measures to, the Licensed Material if doing so
# restricts exercise of the Licensed Rights by any recipient of the Licensed
# Material. No endorsement. Nothing in this Public License constitutes or may be
# construed as permission to assert or imply that You are, or that Your use of the
# Licensed Material is, connected with, or sponsored, endorsed, or granted
# official status by, the Licensor or others designated to receive attribution as
# provided in Section 3(a)(1)(A)(i). Other rights.
#
# Moral rights, such as the right of integrity, are not licensed under this Public
# License, nor are publicity, privacy, and/or other similar personality rights;
# however, to the extent possible, the Licensor waives and/or agrees not to assert
# any such rights held by the Licensor to the limited extent necessary to allow
# You to exercise the Licensed Rights, but not otherwise. Patent and trademark
# rights are not licensed under this Public License. To the extent possible, the
# Licensor waives any right to collect royalties from You for the exercise of the
# Licensed Rights, whether directly or through a collecting society under any
# voluntary or waivable statutory or compulsory licensing scheme. In all other
# cases the Licensor expressly reserves any right to collect such royalties.
# Section 3 – License Conditions.
#
# Your exercise of the Licensed Rights is expressly made subject to the following
# conditions.
#
# Attribution.
#
# If You Share the Licensed Material (including in modified form), You must:
#
# retain the following if it is supplied by the Licensor with the Licensed
# Material: identification of the creator(s) of the Licensed Material and any
# others designated to receive attribution, in any reasonable manner requested by
# the Licensor (including by pseudonym if designated); a copyright notice; a
# notice that refers to this Public License; a notice that refers to the
# disclaimer of warranties; a URI or hyperlink to the Licensed Material to the
# extent reasonably practicable; indicate if You modified the Licensed Material
# and retain an indication of any previous modifications; and indicate the
# Licensed Material is licensed under this Public License, and include the text
# of, or the URI or hyperlink to, this Public License. You may satisfy the
# conditions in Section 3(a)(1) in any reasonable manner based on the medium,
# means, and context in which You Share the Licensed Material. For example, it may
# be reasonable to satisfy the conditions by providing a URI or hyperlink to a
# resource that includes the required information. If requested by the Licensor,
# You must remove any of the information required by Section 3(a)(1)(A) to the
# extent reasonably practicable. If You Share Adapted Material You produce, the
# Adapter's License You apply must not prevent recipients of the Adapted Material
# from complying with this Public License. Section 4 – Sui Generis Database
# Rights.
#
# Where the Licensed Rights include Sui Generis Database Rights that apply to Your
# use of the Licensed Material:
#
#   for the avoidance of doubt, Section 2(a)(1) grants You the right to extract,
# reuse, reproduce, and Share all or a substantial portion of the contents of the
# database; if You include all or a substantial portion of the database contents
# in a database in which You have Sui Generis Database Rights, then the database
# in which You have Sui Generis Database Rights (but not its individual contents)
# is Adapted Material; and You must comply with the conditions in Section 3(a) if
# You Share all or a substantial portion of the contents of the database. For the
# avoidance of doubt, this Section 4 supplements and does not replace Your
# obligations under this Public License where the Licensed Rights include other
# Copyright and Similar Rights. Section 5 – Disclaimer of Warranties and
# Limitation of Liability.
#
# Unless otherwise separately undertaken by the Licensor, to the extent possible,
# the Licensor offers the Licensed Material as-is and as-available, and makes no
# representations or warranties of any kind concerning the Licensed Material,
# whether express, implied, statutory, or other. This includes, without
# limitation, warranties of title, merchantability, fitness for a particular
# purpose, non-infringement, absence of latent or other defects, accuracy, or the
# presence or absence of errors, whether or not known or discoverable. Where
# disclaimers of warranties are not allowed in full or in part, this disclaimer
# may not apply to You. To the extent possible, in no event will the Licensor be
# liable to You on any legal theory (including, without limitation, negligence) or
# otherwise for any direct, special, indirect, incidental, consequential,
# punitive, exemplary, or other losses, costs, expenses, or damages arising out of
# this Public License or use of the Licensed Material, even if the Licensor has
# been advised of the possibility of such losses, costs, expenses, or damages.
# Where a limitation of liability is not allowed in full or in part, this
# limitation may not apply to You. The disclaimer of warranties and limitation of
# liability provided above shall be interpreted in a manner that, to the extent
# possible, most closely approximates an absolute disclaimer and waiver of all
# liability. Section 6 – Term and Termination.
#
# This Public License applies for the term of the Copyright and Similar Rights
# licensed here. However, if You fail to comply with this Public License, then
# Your rights under this Public License terminate automatically. Where Your right
# to use the Licensed Material has terminated under Section 6(a), it reinstates:
#
#   automatically as of the date the violation is cured, provided it is cured within
# 30 days of Your discovery of the violation; or upon express reinstatement by the
# Licensor. For the avoidance of doubt, this Section 6(b) does not affect any
# right the Licensor may have to seek remedies for Your violations of this Public
# License. For the avoidance of doubt, the Licensor may also offer the Licensed
# Material under separate terms or conditions or stop distributing the Licensed
# Material at any time; however, doing so will not terminate this Public License.
# Sections 1, 5, 6, 7, and 8 survive termination of this Public License. Section 7
# – Other Terms and Conditions.
#
# The Licensor shall not be bound by any additional or different terms or
# conditions communicated by You unless expressly agreed. Any arrangements,
# understandings, or agreements regarding the Licensed Material not stated herein
# are separate from and independent of the terms and conditions of this Public
# License. Section 8 – Interpretation.
#
# For the avoidance of doubt, this Public License does not, and shall not be
# interpreted to, reduce, limit, restrict, or impose conditions on any use of the
# Licensed Material that could lawfully be made without permission under this
# Public License. To the extent possible, if any provision of this Public License
# is deemed unenforceable, it shall be automatically reformed to the minimum
# extent necessary to make it enforceable. If the provision cannot be reformed, it
# shall be severed from this Public License without affecting the enforceability
# of the remaining terms and conditions. No term or condition of this Public
# License will be waived and no failure to comply consented to unless expressly
# agreed to by the Licensor. Nothing in this Public License constitutes or may be
# interpreted as a limitation upon, or waiver of, any privileges and immunities
# that apply to the Licensor or You, including from the legal processes of any
# jurisdiction or authority.

# Step 1: download publicly available replication data for Carvalho et al. from https://www.aeaweb.org/articles?id=10.1257/aer.20140481
# Step 2: download publicly available replication code for Farbmacher et al. from http://www.farbmacher.de/codes/FKS_2019_replication.zip
# Step 3: run "KP_replication_script.R" from Farbmacher et al. to produce the data.frame "full.kp.mod"

data = full.kp.mod # See step 3 above

Y1 = data[, 1] # Correct answers per second
Y2 = data[, 2] # Number of correct answers, lower outcome = worse
Y3 = data[, 3] # Total response time, a higher outcome = worse
W = data[, 4] # W = 1 if survey before payday
X = data[, -(1:4)]

# Table 1 Farbmacher
mean(Y2)
sd(Y2)
mean(W)
sd(W)
# Table 2 Farbmacher
print(as.matrix(apply(X, 2, mean)), digits=2)

# Rename variables
colnames(X) = c(
  "live.paycheck.by.paycheck",
  "caloric.crunch",
  "liquidity.constrained",
  "financial.hardship",
  "current.income",
  "pay.day.amount.fraction.of.income",
  "age",
  "male",
  "married",
  "divorced",
  "widowed",
  "never.married",
  "white",
  "black",
  "hispanic",
  "other.race",
  "working",
  "unemployed",
  "disabled",
  "retired",
  "other.employment.status",
  "household.size",
  "college.graduate",
  "some.college",
  "high.school",
  "less.high.school",
  "income_less_5k",
  "income_btw_5k_10k",
  "income_btw_10k_15k",
  "income_btw_15k_20k",
  "income_btw_20k_25k",
  "income_btw_25k_30k",
  "income_btw_30k_35k",
  "income_btw_35k_40k",
  "household.head",
  "children.in.household",
  "metropolitan.area"
)

# Re-encode the income-dummies to an ordinal level (income range) to make it easier for forests to split
inc.dummies = c("income_less_5k", "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k", "income_btw_20k_25k",
                "income_btw_25k_30k", "income_btw_30k_35k", "income_btw_35k_40k")
inc.colnum = which(colnames(X) %in% inc.dummies)
inc.dummies.income = c(5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000)
all(rowSums(X[inc.dummies]) == 1)
income = inc.dummies.income[max.col(X[inc.dummies])]
X = X[-inc.colnum]
X[, "household.income"] = income

# Do the same with education dummies
edu.dummies = c("college.graduate", "some.college", "high.school", "less.high.school")
edu.colnum = which(colnames(X) %in% edu.dummies)
edu.level = c(4, 3, 2, 1)
all(rowSums(X[edu.dummies]) == 1)
edu = edu.level[max.col(X[edu.dummies])] # high = college, low = less than high school
X = X[-edu.colnum]
X[, "education"] = edu

# No need to split on at several digits of precision
X$pay.day.amount.fraction.of.income = round(X$pay.day.amount.fraction.of.income, 2)

# Remove some redundant (perfectly colinear) variables
summary(lm(Y2 ~as.matrix(X)))
X$never.married=NULL
X$other.race=NULL
X$other.employment.status=NULL

out = cbind(outcome.correct.ans.per.second = Y1, outcome.num.correct.ans = Y2, outcome.response.time = Y3,
            treatment = W, X)
write.csv(out, "carvalho2016.csv", row.names = FALSE)
