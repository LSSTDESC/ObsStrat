export TEXINPUTS:=./desc-tex/styles/:./desc-tex/logos/:
export BSTINPUTS:=./desc-tex/bib/:./desc-tex/bst/:

all:
	[ -e .logos ] || { ln -s ./desc-tex/logos .logos >/dev/null 2>&1; }
	latexmk -pdf -g main.tex

clean:
	rm -f *.aux *.fls *.log *Notes.bib *.fdb_latexmk
