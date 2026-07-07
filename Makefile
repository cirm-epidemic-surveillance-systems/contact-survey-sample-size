# Build LaTeX supplements for Overleaf from the Quarto source documents.
#
#   make supplement   Render R/analysis_s1_binning_methods.qmd to a standalone
#                     .tex in overleaf/ -- code and execution are switched off
#                     (the figures already exist), inline math is preserved, and
#                     the figure paths are rewritten to match the Overleaf layout
#                     (figures live in an Overleaf-side supp_figures/ folder). The
#                     figures are copied alongside so overleaf/ compiles as-is.
#
# Workflow: edit the .qmd -> `make supplement` -> push overleaf/ to Overleaf.

QUARTO ?= $(shell command -v quarto 2>/dev/null || echo /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto)

SUP_NAME := analysis_s1_binning_methods
SUP_QMD  := R/$(SUP_NAME).qmd
OUT_DIR  := overleaf
FIGS     := binning_illustration.png binning_moment_convergence.png binning_eigenvalue_convergence.png

.PHONY: supplement
supplement:
	@mkdir -p $(OUT_DIR)/supp_figures
	$(QUARTO) render $(SUP_QMD) --to latex --no-execute -M echo:false
	@# rewrite figure paths (../figures/ -> supp_figures/) for the Overleaf layout
	@sed 's#\.\./figures/#supp_figures/#g' R/$(SUP_NAME).tex > $(OUT_DIR)/$(SUP_NAME).tex
	@rm -f R/$(SUP_NAME).tex
	@cp $(addprefix figures/,$(FIGS)) $(OUT_DIR)/supp_figures/
	@echo "Built $(OUT_DIR)/$(SUP_NAME).tex  (+ supp_figures/)"
