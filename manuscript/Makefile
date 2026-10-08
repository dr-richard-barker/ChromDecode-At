# Build the manuscript and supplement PDFs.
#
#   make          manuscript.pdf + supplementary.pdf (latexmk)
#   make draft    also copies them to ChromDecode-At_draft.pdf / _supplement_draft.pdf
#   make check    every \cite resolves and nothing in references.bib is uncited
#   make clean    remove build artefacts

.PHONY: all draft check clean

all: supplementary.pdf manuscript.pdf

# manuscript.tex reads supplementary.aux (xr-hyper) for Fig. S<n> numbers, so the supplement builds first.

manuscript.pdf: manuscript.tex references.bib supplementary.pdf $(wildcard ../figures/*.png)
	latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript.tex

supplementary.pdf: supplementary.tex $(wildcard ../figures/*.png)
	latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary.tex

draft: all
	cp manuscript.pdf ChromDecode-At_draft.pdf
	cp supplementary.pdf ChromDecode-At_supplement_draft.pdf

check:
	@python3 -c "import re;s=open('manuscript.tex').read();c={k.strip() for g in re.findall(r'\\\\citep?\{([^}]*)\}',s) for k in g.split(',')};b=set(re.findall(r'^@\w+\{([^,]+),',open('references.bib').read(),re.M));print('undefined:',sorted(c-b) or 'none');print('uncited:',sorted(b-c) or 'none');exit(bool(c-b or b-c))"

clean:
	latexmk -C manuscript.tex supplementary.tex 2>/dev/null || true
