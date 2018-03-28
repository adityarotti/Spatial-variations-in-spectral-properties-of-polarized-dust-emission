article.pdf : article.tex #figs/*.pdf
	pdflatex -halt-on-error article.tex
	bibtex article
	pdflatex -halt-on-error article.tex
	pdflatex -halt-on-error article.tex
