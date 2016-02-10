
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{wap}
\LoadClass[10pt,a4paper,onecolumn,oneside]{article}
\usepackage{geometry}
\geometry{left=1.5cm,right=1.5cm,top=2.5cm,bottom=2.5cm}
\usepackage[center]{titlesec}

\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhead{} % clear all fields
\lhead {Journal of Algorithms and Optimization}
\rhead {Jan. 2015, Vol. 3, No. 1, pp. 1-13}

\newcommand*{\TitleFont}{%
      \usefont{\encodingdefault}{cmr}{b}{n}%
      \fontsize{16}{20}%
      \selectfont}
\newcommand*{\secondFont}{%
      \usefont{\encodingdefault}{cmr}{b}{n}%
      \fontsize{14}{16}%
      \selectfont}
\renewcommand\thesection{\Roman{section}}
