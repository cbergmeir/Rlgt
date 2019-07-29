#' @title Getting started with the Rlgt package
#' @description An implementation of LGT and SGT models as 
#' described in ....

#' %%%%%%%%%%%%  Generated using docx2latex.com  %%%%%%%%%%%%%%

#' %%%%%%%%%%%%  v2.0.0-beta  %%%%%%%%%%%%%%

#' \documentclass[12pt]{article}
#' \usepackage{amsmath}
#' \usepackage{latexsym}
#' \usepackage{amsfonts}
#' \usepackage[normalem]{ulem}
#' \usepackage{array}
#' \usepackage{amssymb}
#' \usepackage{graphicx}

#' % If you are compiling in your own LaTeX editor and this 
 #' %part of code is throwing error then remove following few lines of code before \usepackage{subfig}
#' \usepackage[backend=biber,
#' style=numeric,
#' sorting=ynt
#' ]{biblatex}\addbibresource{bibliography.bib}

#' \usepackage{subfig}
#' \usepackage{wrapfig}
#' \usepackage{wasysym}
#' \usepackage{enumitem}
#' \usepackage{adjustbox}
#' \usepackage{ragged2e}
#' \usepackage[svgnames,table]{xcolor}
#' \usepackage{tikz}
#' \usepackage{longtable}
#' \usepackage{changepage}
#' \usepackage{setspace}
#' \usepackage{hhline}
#' \usepackage{multicol}
#' \usepackage{tabto}
#' \usepackage{float}
#' \usepackage{multirow}
#' \usepackage{makecell}
#' \usepackage{fancyhdr}
#' \usepackage[toc,page]{appendix}
#' \usepackage[paperheight=11.69in,paperwidth=8.26in,left=1.25in,right=1.25in,top=1.0in,bottom=1.0in,headheight=1in]{geometry}
#' \usepackage[utf8]{inputenc}
#' \usepackage[T1]{fontenc}
#' \usepackage[hidelinks]{hyperref}
#' \usetikzlibrary{shapes.symbols,shapes.geometric,shadows,arrows.meta}
#' \tikzset{>={Latex[width=1.5mm,length=2mm]}}
#' \usepackage{flowchart}\TabPositions{0.5in,1.0in,1.5in,2.0in,2.5in,3.0in,3.5in,4.0in,4.5in,5.0in,5.5in,}

#' \urlstyle{same}


 #' %%%%%%%%%%%%  Set Depths for Sections  %%%%%%%%%%%%%%

#' % 1) Section
#' % 1.1) SubSection
#' % 1.1.1) SubSubSection
#' % 1.1.1.1) Paragraph
#' % 1.1.1.1.1) Subparagraph


#' \setcounter{tocdepth}{5}
#' \setcounter{secnumdepth}{5}


 #' %%%%%%%%%%%%  Set Depths for Nested Lists created by \begin{enumerate}  %%%%%%%%%%%%%%


#' \setlistdepth{9}
#' \renewlist{enumerate}{enumerate}{9}
	#' \setlist[enumerate,1]{label=\arabic*)}
	#' \setlist[enumerate,2]{label=\alph*)}
	#' \setlist[enumerate,3]{label=(\roman*)}
	#' \setlist[enumerate,4]{label=(\arabic*)}
	#' \setlist[enumerate,5]{label=(\Alph*)}
	#' \setlist[enumerate,6]{label=(\Roman*)}
	#' \setlist[enumerate,7]{label=\arabic*}
	#' \setlist[enumerate,8]{label=\alph*}
	#' \setlist[enumerate,9]{label=\roman*}

#' \renewlist{itemize}{itemize}{9}
	#' \setlist[itemize]{label=$\cdot$}
	#' \setlist[itemize,1]{label=\textbullet}
	#' \setlist[itemize,2]{label=$\circ$}
	#' \setlist[itemize,3]{label=$\ast$}
	#' \setlist[itemize,4]{label=$\dagger$}
	#' \setlist[itemize,5]{label=$\triangleright$}
	#' \setlist[itemize,6]{label=$\bigstar$}
	#' \setlist[itemize,7]{label=$\blacklozenge$}
	#' \setlist[itemize,8]{label=$\prime$}



 #' %%%%%%%%%%%%  Header here  %%%%%%%%%%%%%%


#' \pagestyle{fancy}
#' \fancyhf{}
#' \cfoot{ 
#' \vspace{\baselineskip}

#' \vspace{\baselineskip}
#' }
#' \renewcommand{\headrulewidth}{0pt}
#' \setlength{\topsep}{0pt}\setlength{\parskip}{9.96pt}
#' \setlength{\parindent}{0pt}
#' \renewcommand{\arraystretch}{1.3}


#' %%%%%%%%%%%%%%%%%%%% Document code starts here %%%%%%%%%%%%%%%%%%%%



#' \begin{document}
#' \subsubsection*{4.3.1. LGT (Local and Global Trend)}
#' \addcontentsline{toc}{subsubsection}{4.3.1. LGT (Local and Global Trend)}
#' In terms of mathematical notation, the model can be fully represented as follow:\par

 #' \( y_{t+1} \sim Student \left( v,y_{t+1}, \sigma _{t+1} \right)  \)  \  \ \ \ \ \ \ \ \ \  \   \ \  (eq.1.1)\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \par

 #' \( y_{t+1}=l_{t}+ \gamma l_{t}^{ \rho }+ \lambda b_{t} \) \ \ \ \  \tab \tab \ \ \ \ \  (eq. 1.2)\tab \tab \tab \par

 #' \( l_{t}= \alpha y_{t}+ \left( 1- \alpha  \right)  \left(  l_{t-1} \right)  \) \tab \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  (eq. 1.3)\tab \par

 #' \( b_{t+1}= \beta  \left( l_{t+1}-l_{t} \right) + \left( 1- \beta  \right) b_{t} \) \ \ \ \ \ \ \ \ \ \ \ \ \  (eq. 1.4)\tab \tab \par

 #' \(  \sigma _{t+1}= \sigma l_{t}^{ \tau}+ \varsigma  \) \tab \tab \tab \ \ \ \ \ \ \  (eq. 1.5)\tab \tab \tab \par

#' and the notations are defined as follow:\par

#' \uline{Time-varying variables}\par

#' \setlength{\parskip}{0.0pt}
 #' \( y_{t} \) \tab : value of the dependent variable of interest at time t\par

 #' \( y_{t+1} \) \tab : predicted value of y at time t+1 given information up to time t\par

 #' \(  \sigma _{t} \)  \tab : variance of the distribution at time t\par

 #' \( l_{t} \) \tab : level at time t\par

#' \setlength{\parskip}{6.0pt}
 #' \( b_{t} \) \tab : local trend at time t\par

#' \uline{Parameters of the model which need to be estimated}\par

#' \setlength{\parskip}{0.0pt}
 #' \( v \) \tab : degrees of freedom of the t-distribution\par

#' \setlength{\parskip}{6.0pt}
 #' \(  \gamma  \) \tab : coefficient of the global trend\par

 #' \(  \rho  \) \tab : power coefficient of the global trend\par

 #' \(  \lambda  \) \tab : damping coefficient of the local trend\par

 #' \(  \alpha  \) \tab : smoothing parameter for the level\par

 #' \(  \beta  \) \tab : smoothing parameter for the local trend\par

 #' \(  \sigma  \) \tab : coefficient of the heteroscedastic size of error function\par

 #' \(  \tau \) \tab : power coefficient of the heteroscedastic size of error function\par

 #' \(  \varsigma  \) \tab : base/ minimum value of the size of error function\par

#' \setlength{\parskip}{9.96pt}
#' The rationale of each individual equation of the model is discussed below.\par

#' \chapter{eq. 1.1. Student’s t error distribution}\par

#' The data value follows Student’s t-distribution around the expected value of the data with a time-varying standard deviation. The Student-t distribution can be seen as a generalisation of the normal distribution to allow for a fat-tailed error distribution \par

#' \chapter{eq. 1.5. Heteroscedastic Error }\par

#' In addition to accounting for possible fat-tailed error distribution, the error function\ also allows the variance of the error to change as the level changes.  This is achieved by allowing the scale (deviation) parameter of the assumed Student’s t-distribution to vary in proportion to the current level of the time series. This will account for common situations, where the magnitude of the error will increase as the value of the data points increases. \par

#' Moreover, this relationship does not necessarily be linear as the parameter  \(  \tau \)  controls the growth of the variance of the error term. In practice, the values taken by the parameter is often limited between 0 and 1. A value of 0 corresponds to constant variance, i.e. homoscedasticity, whereas a value of 1 approximates the behaviour of the multiplicative error ETS model. The value between 0 and 1 will describe an error function which grows in relation to the increase in data value, but at a slower pace than linear growth.\par

#' \chapter{eq. 1.2: One-step prediction forecast}\par

#' \setlength{\parskip}{5.04pt}
#' There are three distinct terms that constitute this Bayesian ETS model: a level term, and a couple of different trends. The first term \(  l_{t}  \) is the level term, while the second term  \( l_{t}^{ \rho } \)  refers to the global trend which increases with the level of the dependent variables  \( l_{t} \)  at a constant $``$rate$"$   \(  \gamma  \) . Similar to the heteroscedasticity in eq. 1.5, the change in the value of the dependant variable does not necessarily grow linearly with respect to the level of  \( the variable. \)  This relationship can be tuned in by the use of exponential parameter  \(  \rho  \) . The interpretation of this global trend is also analogous to the time-varying error term, i.e. the value of  \( 0< \rho <1 \)  indicates a global trend which grows faster than the additive models but slower than the multiplicative model.\par

#' The last term of the right-hand side of eq.1.2. refers to the usual local dampen trend in ETS model. However, there is an addition of dampen parameter  \(  \lambda  \) , constrained such that  \( -1<  \lambda  <1 \) , to reduce the strength of the local trend model. \par

#' \setlength{\parskip}{9.96pt}
#' \chapter{eq. 1.3: Level adjustment equation}\par

#' \setlength{\parskip}{5.04pt}
#' This equation is defined according to the classical linear trend ETS model. The level at time t ( \( l_{t} \) ) is calculated as a weighted average of the current observation  \( y_{t} \)  and the previous level at lag 1 ( \( l_{t-1} \) ) with smoothing parameter  \(  \alpha  \) .\par

#' \setlength{\parskip}{9.96pt}
#' \chapter{Eq. 1.4: Local trend adjustment equation}\par

#' \setlength{\parskip}{5.04pt}
#' Similarly, the evolution of the local trend  \( b_{t}  \) is identically defined to the linear trend method.\  The local trend at time t  \( b_{t} \)  is obtained as a weighted average of the difference in level terms  \(  \left( l_{t}-l_{t-1} \right)  \)  and the trend at time t-1 ( \( b \_  \left( t-1 \right)  \right)  \)  with the smoothing parameter  \(  \beta  \) .\par

#' \setlength{\parskip}{9.96pt}
#' \subsubsection*{4.3.3. Parameters’ Prior Distributions in L/SGT Models}
#' \addcontentsline{toc}{subsubsection}{4.3.3. Parameters’ Prior Distributions in L/SGT Models}
#' Although the exact prior distributions for the parameters of the L/SGT models are not exactly specified in the technical report, it is likely that Smyl used the following prior distributions based on his original code in Rlgt package (refer to chapter 5):\par

#' Note that a more elaborate explanation of these prior choices will also be given in Chapter 5.\par

#' \setlength{\parskip}{6.0pt}
#' \begin{adjustwidth}{0.69in}{0.0in}
 #' \(  \sigma  \) ,  \(  \gamma ,  \varsigma  \) \tab : Cauchy distribution with 0 location value and the scale parameter equals to 1/200 of the maximum value of y\par

#' \end{adjustwidth}

#' \begin{adjustwidth}{0.69in}{0.0in}
#' b\tab : Normally distributed with a mean of 0 and standard deviation of 1/200 of the maximum value of y.\par

#' \end{adjustwidth}

#' \begin{adjustwidth}{0.69in}{0.0in}
 #' \(  \phi  \) \  \tab : Uniformly distributed between -1 and 1.\par

#' \end{adjustwidth}

#' \begin{adjustwidth}{0.69in}{0.0in}
 #' \(  \alpha  \) , \(  \beta  \) ,  \(  \zeta  \) \tab : Uniform between 0 and 1.\par

#' \end{adjustwidth}

#' \begin{adjustwidth}{0.69in}{0.0in}
 #' \( s_{t} \) \tab : (i.i.d) normal with a mean of 1 and standard deviation of 0.3 before being normalised.\par

#' \end{adjustwidth}

#' \setlength{\parskip}{9.96pt}
#' \begin{adjustwidth}{0.69in}{0.0in}
 #' \(  \rho  \) \tab : Uniform between -0.5 and 1.0\par

#' \end{adjustwidth}

#' \begin{adjustwidth}{0.69in}{0.0in}
 #' \( v \) \tab : Uniformly distributed between 2 and 20\par

#' \end{adjustwidth}

#' \begin{adjustwidth}{0.69in}{0.0in}
 #' \(  \tau \)  \tab : Beta distribution with shape parameters  \(  \alpha =1 \)  and \(   \beta =1 \) \par

#' \end{adjustwidth}


#' \vspace{\baselineskip}

#' \vspace{\baselineskip}

#' \vspace{\baselineskip}
#' \end{document}

#' @name Rlgt-package
#' @aliases Rlgt
#' @docType package
#' @title Getting started with the Rlgt package
#' @keywords forecasting, exponential smoothing
# @exportPattern "^[[:alpha:]]+"
#' @useDynLib Rlgt, .registration=TRUE
#' @import methods
#' @import Rcpp
#'
#' @examples
#' x <- 1
#' @docType package
#' @name Rlgt
NULL

