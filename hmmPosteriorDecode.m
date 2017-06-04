(* decode
INPUT
 - observationSeq is a list of observations, e.g., {2,3,4,1,2}
 - states is a list of the state names, e.g., {m, h}
 - alphabet is a list of the HMM alphabet, e.g., {1, 2, 3, 4}
 - emissionMatrix is a matrix of dimensions {Length[states], Length[alphabet]}.  
     emissionMatrix[[i,j]] is the probability of emitting letter j from state i, 
     e.g., {{0.4, 0.1, 0.1, 0.4}, {0.05, 0.4, 0.5, 0.05}}
 - transitionMatrix is a matrix of dimensions {Length[states], Length[states]}.
     transitionMatrix[[i, j]] is the probability of transitioning to state j on one transition starting from state i.
     e.g., {{0.99, 0.01}, {0.01, 0.99}}
 - initialStateProbs is a list of dimensions {Length[states]}
     initialStateProbs[[i]] is the prior probability that the state from which the first observations was
     emitted is state i.  
OUTPUT
- stateSeq is a list of dimensions {Length[observationSeq]}.
  stateSeq[[i]] is the ith state in the the most likely sequence of states, given the observations. 
  e.g., {h,h,m,m,m}.
  *)

(* TODO: Remove the comments surround this decode function template, complete the function, and test it. *)
decode[observationSeq_, {states_, alphabet_, emissionMatrix_, transitionMatrix_, initialStateProbs_}] := 
	Module[{stateSeq={},maxViterbiTable={}, viterbiTable ={},tempViterbiList={},tempList={},k,i,j,m,l,p,s,q},
		
		(*add the maxium for the first column*)
		Do[AppendTo[tempViterbiList, initialStateProbs[[k]]*emissionMatrix[[k, observationSeq[[1]]]]],{k,1,Length[emissionMatrix]}];
		AppendTo[viterbiTable,Split[tempViterbiList]];
		tempViterbiList={};
		
		(*calculate the rest of the maxium based on the oberseration*)
		For[i=2,i<=Length[observationSeq],i++,
			tempViterbiList={};
			For[j=1,j<=Length[states],j++,
				Do[AppendTo[tempList, transitionMatrix[[m,j]]*emissionMatrix[[j,observationSeq[[i]]]]*Max[viterbiTable[[i-1,m]]]
					], {m,1,Length[states]}];
				AppendTo[tempViterbiList, tempList];
				tempList={};
			];
			AppendTo[viterbiTable,tempViterbiList];
		];
		
		(*take the maxium for each column - observation *)
		Do[AppendTo[maxViterbiTable,Max[viterbiTable[[s]]] ],{s,1,Length[observationSeq]}];
		
		For[p=1,p<=Length[states],p++,
					If[viterbiTable[[1,p,1]]==maxViterbiTable[[1]],AppendTo[stateSeq,states[[p]]]]
		];
		
		(*based on the maxium, find the most likely state sequence*)
		For[l=2,l<=Length[observationSeq],l++,
			For[p=1,p<=Length[states],p++,
				For[q=1,q<=Length[states],q++,
					If[viterbiTable[[l,p,q]]==maxViterbiTable[[l]],AppendTo[stateSeq,states[[p]]]]
				]
				
			]
		];
		Return[stateSeq];
		
		(* Return the state sequence *)
	];
	

(*********************************
*posterior decoding
***********************************)
	
(*function to calculate alphas with forward algorithm*)
forward[observationSeq_, {states_, alphabet_, emissionMatrix_, transitionMatrix_, initialStateProbs_}] := 
	Module[{viterbiTable ={},tempViterbiList={},tempList={},k,i,j,m},
		
		(*add the maxium for the first column *)
		Do[AppendTo[tempViterbiList, initialStateProbs[[k]]*emissionMatrix[[k, observationSeq[[1]]]]],{k,1,Length[emissionMatrix]}];
		AppendTo[viterbiTable,tempViterbiList];
		tempViterbiList={};
		
		(*calculate the rest of the alphas based on the oberseration*)
		For[i=2,i<=Length[observationSeq],i++,
			tempViterbiList={};
			For[j=1,j<=Length[states],j++,
				Do[AppendTo[tempList, transitionMatrix[[m,j]]*viterbiTable[[i-1,m]]
					], {m,1,Length[states]}];
				AppendTo[tempViterbiList, emissionMatrix[[j,observationSeq[[i]]]]*Total[tempList]];
				tempList={};
			];
			AppendTo[viterbiTable,tempViterbiList];
		];	
		Return[viterbiTable];	
	];

(*function to calculate betas with backward algorithm*)
backward[observationSeq_, {states_, alphabet_, emissionMatrix_, transitionMatrix_, initialStateProbs_}] :=
	Module[{tempBetaList ={},betaList={},tempList={},i,j,k,m },
		(*add 1 to beta List for the last observation*)
		Do[AppendTo[tempBetaList, 1],{k,1,Length[states]}];
		AppendTo[betaList,tempBetaList];
		tempBetaList={};
		
		(*calculate the rest of betas*)
		For[i=Length[observationSeq]-1,i>=1,i--,
			tempBetaList={};
			For[j=1, j<=Length[states],j++,
				Do[AppendTo[tempList,transitionMatrix[[j,m]]*emissionMatrix[[m,observationSeq[[i+1]]]]*betaList[[Length[observationSeq]-i,m]]],
					{m,1,Length[states]}];
				AppendTo[tempBetaList, Total[tempList]];
				tempList={};
			];
			AppendTo[betaList,tempBetaList];
		];
		Return[Reverse[betaList]];
	]
(*function to calculate the posterior based on the data from forward and backward algorithm*)
posterior[forwardList_,backwardList_]:=
	Module[{postList={},subSum={},posteriors={},i},
		postList=forwardList*backwardList;
		Do[AppendTo[subSum,Total[postList[[i]]]],{i,1,Length[postList]}];
		posteriors = postList/subSum;
		Return[posteriors];
	]
(*decode the sequence states based on the posterior probability of state*)
posteriorDecode[observationSeq_, {states_, alphabet_, emissionMatrix_, transitionMatrix_, initialStateProbs_}] := 
	Module[{forwardList={},backwardList={},posteriorList={},stateSeq={},l,p},
		forwardList=forward[observationSeq, {states, alphabet, emissionMatrix, transitionMatrix, initialStateProbs}];
		backwardList =backward[observationSeq, {states, alphabet, emissionMatrix, transitionMatrix, initialStateProbs}];
		posteriorList =posterior[forwardList,backwardList];
		
		(*take the maxium for each row - observation *)
		(*Do[AppendTo[maxPosterior,Max[posteriorList[[i]]] ],{i,1,Length[observationSeq]}];*)
		
		(*based on the maxium, find the most likely state sequence*)
		For[l=1,l<=Length[observationSeq],l++,
			For[p=1,p<=Length[states],p++,
				If[posteriorList[[l,p]]==Max[posteriorList[[l]]],AppendTo[stateSeq,states[[p]]]]
			];
		];
		Return[stateSeq];
	]
(* calculateAccuracy takes a state sequence genereted from mixed2.fa and calculates 
the number of correctly labeled states.  Note: this function only computes the
accuracy for the mixed2.fa observations.

INPUT
 - stateSeq is a list of state sequences, e.g., {h,m,h,m,m}
 
 OUTPUT
 - numCorrectStates = [int] number of correcly labeled states.

*)
	
calculateAccuracy[stateSeq_] := 
	Module[{keyStateSequence, numCorrectStates},
	
	keyStateSequence = Flatten[Characters[ToLowerCase[Import["mixed2key.fa"]]]];
	numCorrectStates = Count[MapThread[Equal, {stateSeq, keyStateSequence}], True]
	];

(* readHMM takes a HMM text file and outputs the state names, the alphabet
the transition matrix, and emission matrix of the HMM

INPUT
 - file is a path to an HMM file (see notebook for the format of HMM files). 

  OUTPUT
 - states = [list] list of the state names, e.g., {m, h}
 - alphabet = [list] list of the HMM alphabet, e.g., {1, 2, 3, 4}
 - emissionMatrix = [matrix of size numStates x numAlphabet] the emission matrix.  
     Element eij = the probability of state i emitting letter j., e.g., {{0.4, 0.1, 0.1, 0.4}, {0.05, 0.4, 0.5, 0.05}}
 - transitionMatrix = [matrix of size numStates x numStates] the transition matrix.
     Element tij = the probability of state i transitioning to state j, e.g., {{0.99, 0.01}, {0.01, 0.99}}

*)

(*Note: this is not exactly how I would hav written readHMM stylewise, but it works so I won't change it for now. -MRB *)

readHMM[file_] := 
	Module[{a, numStates, alphabet, numAlphabet, firstStateIndex, lastStateIndex,
		states, firstStateProbIndex, lastStateProbIndex, initialStateProbs, 
		firstEmissionIndex, lastEmissionIndex, emissionList, emissionMatrix,
		firstTransitionIndex, lastTransitionIndex, transitionList, transitionMatrix}, 
		
	a = Import[file, {"Text", "Words"}];
	
	numStates = ToExpression[a[[1]]]; (* Use ToExpression to convert from character to number *)

	alphabet = Characters[a[[2]]];
	numAlphabet = Length[alphabet];

	firstStateIndex = 3;
	lastStateIndex = firstStateIndex + numStates - 1;
	states = a[[firstStateIndex ;; lastStateIndex]];

	firstStateProbIndex = lastStateIndex + 1;
	lastStateProbIndex = firstStateProbIndex + numStates - 1;
	initialStateProbs = ToExpression[a[[firstStateProbIndex ;; lastStateProbIndex]]];

	firstEmissionIndex = lastStateProbIndex + 1;
	lastEmissionIndex = firstEmissionIndex + numStates*numAlphabet - 1;
	emissionList = ToExpression[a[[firstEmissionIndex ;; lastEmissionIndex]]];
	emissionMatrix = Partition[emissionList, numAlphabet];

	firstTransitionIndex = lastEmissionIndex + 1;
	lastTransitionIndex = firstTransitionIndex + numStates*numStates - 1;
	transitionList = ToExpression[a[[firstTransitionIndex ;; lastTransitionIndex]]];
	transitionMatrix = Partition[transitionList, numStates];
	
	{states, alphabet, emissionMatrix, transitionMatrix, initialStateProbs}

];

	
(* readFasta reads a fasta file and outputs the nucleotide sequence converted to numbers
INPUT
- fastaFile is a string representing the path to fasta file

OUTPUT
- input is a list of bases in the file indicated by fastaFile.  
  bases are translated form ACGT to 1234.
  e.g., {1,3,2,4,2}
*)
readFasta[fastaFile_]:=
	Flatten[Map[Characters, Import[fastaFile]] 
		   /. {"A"->1, "C"->2, "G"->3, "T"->4}
		   ]