\name{internals}

%%%
%%% Model fitting
%%%
\alias{rcoxFit}
\alias{rcoxFitPrimitive}
\alias{rcoxFitPrimitive.rcon}
\alias{rcoxFitPrimitive.rcor}

\alias{rconFit}
\alias{rcorFit}

\alias{rconFitIterative}
\alias{rcorFitIterative}
\alias{rconFitHYD}
\alias{rcorFitHYD}
\alias{rconFitHYS}

\alias{fitIPSedge}
\alias{fitIPSset}
\alias{fitNR}
\alias{fitHYS}
\alias{modifiedIPS}
\alias{modifiedNewton}



%%%
%%% Incidence matrices etc.
%%%
\alias{makeIncMatList}
\alias{makeIncMat}
\alias{makeIncMatECC}
\alias{makeIncMatVCC}


%%%
%%% Changing model representations between lists and formulae
%%% (Silly names, change these!)
%%%
\alias{joinForm}
\alias{splitForm}

%%%
%%% rcox settings and control
%%%
\alias{setControl}
\alias{summary.rcox}
\alias{findKinModelPrimitive}
\alias{findKinModel.rcon}
\alias{findKinModel.rcor}
\alias{findKinModel}
\alias{copyInitial}
%%%
%%% Numerical functions
%%%

\alias{cholSolve}
\alias{ellK}


%%%
%%% Model handling
%%%

\alias{getCC}
\alias{getIdxOfCC}
\alias{namesToPairs}

%%%
%%% Score, information and scoring algorithm
%%%

\alias{getScore}
\alias{getScore.rcon}
\alias{getScore.rcor}
\alias{getInf}
\alias{getInf.rcon}
\alias{getInf.rcor}
\alias{getKFromTheta}
\alias{getKFromTheta.rcon}
\alias{getKFromTheta.rcor}
\alias{getThetaFromK}
\alias{getThetaFromK.rcon}
\alias{getThetaFromK.rcor}
\alias{multivariateNewton}



%%%
%%% To be sorted
%%%

\alias{cclist2df}
\alias{foooo}
\alias{getParmVecEntry}
\alias{print.cctest}
\alias{printIfDetails}
\alias{showIfDisplay}
\alias{testSplitAllCC}



\alias{cliqueToEdges}

\alias{getCovConc}
\alias{getCovCorr}
\alias{getdf}
\alias{getFormula}
\alias{getIndex}
\alias{getIndex.character}
\alias{getIndex.integer}
\alias{getIndex.list}
\alias{getIndex.numeric}
\alias{getNames}
\alias{getNames.list}
\alias{getNames.numeric}
\alias{getSlot}

\alias{listOrder}
\alias{listOrder.default}
\alias{listOrder.list}
\alias{listOrder.numeric}


\alias{mparse}
\alias{mparsemodel}
\alias{print.rcox}
\alias{rcorfit2}

\alias{regulariseS}
\alias{rescaleC}
\alias{shprint}
 
\alias{rcox.control}

\alias{varDiffConc}
\alias{varDiffCorr}
\alias{str2formula}

\alias{findCanonicalRepresentation}
\alias{findFitRepresentations}
\alias{getAtomCC}
\alias{getCompCC}
\alias{testAddEdges}


\title{Internal functions for the gRcox package}
\description{
  These functions should not be called directly by the user.
}

\author{Søren Højsgaard, sorenh@agrsci.dk}
\keyword{utilities}% at least one, from doc/KEYWORDS

