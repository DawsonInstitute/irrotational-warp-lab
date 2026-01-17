# Paper Figure Generation Pipeline
# Regenerates all figures referenced in docs/paper/main.tex

PAPER_DIR := docs/paper
FIGURES_DIR := $(PAPER_DIR)/figures
RESULTS_DIR := results

# Python environment
PYTHON := .venv/bin/python

.PHONY: all figures paper clean help

all: figures paper

help:
	@echo "Irrotational Warp Lab - Paper Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all       - Build figures and compile paper (default)"
	@echo "  figures   - Regenerate all paper figures"
	@echo "  paper     - Compile LaTeX document"
	@echo "  clean     - Remove generated files"
	@echo "  test      - Run full test suite"
	@echo ""
	@echo "Requirements:"
	@echo "  - Python virtual environment at .venv/"
	@echo "  - LaTeX distribution (pdflatex, bibtex)"
	@echo ""

# Create figures directory
$(FIGURES_DIR):
	mkdir -p $(FIGURES_DIR)

# Figure 1: 3D Convergence Study
$(FIGURES_DIR)/convergence_3d.pdf: $(FIGURES_DIR) scripts/make_paper_figures.py
	@echo "Generating convergence study figure..."
	$(PYTHON) scripts/make_paper_figures.py --figure convergence_3d \
		--out $(FIGURES_DIR)/convergence_3d.pdf

# Figure 2: Superluminal Sweep
$(FIGURES_DIR)/superluminal_sweep.pdf: $(FIGURES_DIR) scripts/make_paper_figures.py
	@echo "Generating superluminal sweep figure..."
	$(PYTHON) scripts/make_paper_figures.py --figure superluminal \
		--out $(FIGURES_DIR)/superluminal_sweep.pdf

# Figure 3: Optimization Comparison
$(FIGURES_DIR)/optimization_comparison.pdf: $(FIGURES_DIR) scripts/make_paper_figures.py
	@echo "Generating optimization comparison figure..."
	$(PYTHON) scripts/make_paper_figures.py --figure optimization \
		--out $(FIGURES_DIR)/optimization_comparison.pdf

# All figures
figures: $(FIGURES_DIR)/convergence_3d.pdf \
         $(FIGURES_DIR)/superluminal_sweep.pdf \
         $(FIGURES_DIR)/optimization_comparison.pdf
	@echo "✓ All figures generated"

# Compile paper
paper: $(PAPER_DIR)/main.pdf

$(PAPER_DIR)/main.pdf: $(PAPER_DIR)/main.tex figures
	@echo "Compiling LaTeX document..."
	cd $(PAPER_DIR) && pdflatex -interaction=nonstopmode main.tex
	cd $(PAPER_DIR) && bibtex main || true
	cd $(PAPER_DIR) && pdflatex -interaction=nonstopmode main.tex
	cd $(PAPER_DIR) && pdflatex -interaction=nonstopmode main.tex
	@echo "✓ Paper compiled: $(PAPER_DIR)/main.pdf"

# Run tests
test:
	$(PYTHON) -m pytest -q

# Clean build artifacts
clean:
	rm -rf $(FIGURES_DIR)/*.pdf
	cd $(PAPER_DIR) && rm -f *.aux *.log *.out *.bbl *.blg *.toc *.pdf
	@echo "✓ Cleaned build artifacts"

# Quick draft (single LaTeX pass, no bibtex)
draft: figures
	cd $(PAPER_DIR) && pdflatex -interaction=nonstopmode main.tex
	@echo "✓ Draft compiled (single pass)"
