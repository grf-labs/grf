Sys.setenv(LINTR_COMMENT_BOT='false')
library(lintr)

linters = with_defaults(
  line_length_linter(120), # Max line length = 120
  object_name_linter = NULL, # To allow variable.names + function_names
  commented_code_linter = NULL, # Misc. false positives
  object_usage_linter = NULL, # Misc. false positives
  spaces_left_parentheses_linter = NULL, # Misc. false positives
  default = default_linters);

res = lint_package(path = 'grf/', linters = linters, exclusions = 'grf/R/RcppExports.R')
res

if (length(res) > 0) quit(status = 1)
