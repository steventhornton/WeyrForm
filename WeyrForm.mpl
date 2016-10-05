# ======================================================================= #
# ======================================================================= #
#                                                                         #
# WeyrForm.mpl                                                            #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 5/2016                                                 #
#                                                                         #
# A function for computing the Weyr canonical form of a matrix            #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   WeyrForm(A)                                                           #
#   WeyrForm(A, output)                                                   #
#                                                                         #
# INPUT                                                                   #
#   A ........ Square matrix                                              #
#   output ... Optional, default: output = W                              #
#              Specify what to return:                                    #
#                   - W:      Return the Weyr form                        #
#                   - Q:      Return the transformation matrix            #
#                   - [W, Q]: Return the Weyr form and transformation     #
#                             matrix                                      #
#                   - [Q, W]: Same as previous just different order.      #
#                                                                         #
# OUTPUT                                                                  #
#                                                                         #
# REFERENCES                                                              #
#   O'Meara, K. (2011). Advanced Topics in Linear Algebra: Weaving Matrix #
#   Problems through the Weyr Form: Weaving Matrix Problems through the   #
#   Weyr Form. Oxford University Press, USA.                              #
# ======================================================================= #
# ======================================================================= #
WeyrForm := module()

    export ModuleApply;

    local
        processInput,
        implementation,
        implementation_W,
        implementation_WQ,
        getJordanStructure,
        eigsAndBlockSizeFromJCF,
        jordanStructureToWeyrStructure,
        weyrBlockMatrix,
        sortJordanForm,
        sortJordanBlock,
        JCF_to_WCF_Transformation_One_Eig,
        permutationMatrix;
        
    
    uses LinearAlgebra, ListTools;

    ModuleApply := proc() :: {Matrix, seq(Matrix, Matrix)};
        local A, out;
        A, out := processInput(args);
        return implementation(A, out);
    end proc;

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# METHODS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# ---------------------------------------------------------------------------- #
# processInput
#
# Cleanse the input values and throw errors if the input values do not take the correct form.
#
# INPUT
#
# OUTPUT
#   
# ---------------------------------------------------------------------------- #
processInput := proc(A::Matrix, out::list(string), $)
    
    local m, n;

    # Check that A is a matrix
    # Check that A is a square matrix
    # Check that A contains Gaussian integers

    # Check that the output is of the form W, Q, [W], [Q], [W, Q] or [Q, W]

    m, n := op(1,A);
    if m <> n then
        error "expecting a square Matrix, got dimensions %1, %2", m, n
    end if;
    
    return A, "W";

end proc:

implementation := proc(A::Matrix, output, $)
    return implementation_W(A);
end proc:


# ----------------------------------------------------------------------- #
# implementation_W                                                        #
#                                                                         #
# Compute the Weyr canonical form of a matrix.                            #
#                                                                         #
# INPUT                                                                   #
#   A ... A square matrix                                                 #
#                                                                         #
# OUTPUT                                                                  #
#   The Weyr form of A                                                    #
# ----------------------------------------------------------------------- #
implementation_W := proc(A::Matrix, $)::Matrix;
    
    local J::Matrix,
        jordanStructure::table,
        weyrStructure::table,
        eigVal,
        weyrBlockList::list(Matrix),
        W::Matrix;
    
    # Compute Jordan form of A
    J := JordanForm(A);
    
    # Get the Weyr structure
    jordanStructure := getJordanStructure(J);
    weyrStructure := map(jordanStructureToWeyrStructure, jordanStructure);
    
    # Generate a list of the Weyr block matrices
    weyrBlockList := [];
    for eigVal in indices(weyrStructure, 'nolist') do
        W := weyrBlockMatrix(eigVal, weyrStructure[eigVal]);
        weyrBlockList := [op(weyrBlockList), W];
    end do;
    
    W := DiagonalMatrix(weyrBlockList);
    
    ASSERT(IsSimilar(W, J), "W must be similar to J");
    ASSERT(IsSimilar(W, A), "W must be similar to A");
    
    return W;
    
end proc;


# ----------------------------------------------------------------------- #
# implementation_WQ                                                       #
#                                                                         #
# Compute the Weyr canonical form and similarity transformation matrix.   #
#                                                                         #
# INPUT                                                                   #
#   A ... A square matrix                                                 #
#                                                                         #
# OUTPUT                                                                  #
#   A sequence of two matrices, the first is the Weyr canonical form of   #
#   the input matrix, and the second is the similarity transformation     #
#   matrix.                                                               #
# ----------------------------------------------------------------------- #
# Return both W and P s.t. W = P^-1 A P
implementation_WQ := proc(A::Matrix, $)

    local J::Matrix, 
          Q1::Matrix,
          Q2::Matrix,
          QQ::Matrix,blockList, 
        startBlock::posint, 
          currentBlockEig, 
          i::posint, 
          permutations::list(Matrix), 
          W::Matrix;
    
    # Compute Jordan form of A
    Q1 := JordanForm(A, 'output' = 'Q');
    
    J := MatrixInverse(Q1).A.Q1;
    
    Q2 := sortJordanForm(J);
    
    # Update Q and J
    QQ := Q1.Q2;
    J := MatrixInverse(QQ).A.QQ;
    
    # Split J by eigenvalue
    # For each Jordan block matrix, compute the transformation to Weyr form
    blockList := [];
    startBlock := 1;
    currentBlockEig := J[1,1];
    for i to RowDimension(A) do
        if currentBlockEig <> J[i,i] then
            blockList := [op(blockList), J[startBlock..(i-1), startBlock..(i-1)]];
            currentBlockEig := J[i,i];
            startBlock := i;
        elif i = RowDimension(A) then
            blockList := [op(blockList), J[startBlock..i, startBlock..i]];
        end if;
    end do;
    
    permutations := map(JCF_to_WCF_Transformation_One_Eig, blockList);
    
    Q1 := DiagonalMatrix(permutations);
    Q1 := MatrixInverse(Q1);
    
    QQ := QQ.Q1;
    
    W := MatrixInverse(QQ).A.QQ;
    
    ASSERT(IsSimilar(A, W), "W and A must be similar");
    
    return W, QQ;
    
end proc;


# Return P s.t. W = P^-1 A P
#implementation_P := proc(J::Matrix, Q::Matrix)
#end proc:


# ----------------------------------------------------------------------- #
# getJordanStructure                                                      #
#                                                                         #
# Determines the Jordan structure of a matrix. That is, given a Jordan    #
# block matrix, it returns a table where the indices are the unique       #
# eigenvalues from the input matrix, and the entries in the table are     #
# lists of positive integers in decreasing order corresponding to the     #
# Jordan blocks in the input matrix for that eigenvalue.                  #
#                                                                         #
# INPUT                                                                   #
#   J ... Jordan block matrix                                             #
#                                                                         #
# OUTPUT                                                                  #
#   A table as described above.                                           #
# ----------------------------------------------------------------------- #
getJordanStructure := proc(J::Matrix, $)::table;
    
    local eigsAndSize::list(list),
                 pair::[anything,posint],
      jordanStructure::table,
                 item;
    
    # Case where J is a 1x1 matrix
    if RowDimension(J) = 1 then
        return table([J[1,1] = 1]);
    end if;
    
    eigsAndSize := eigsAndBlockSizeFromJCF(J);
    
    jordanStructure := table(); 
    
    # Make a table where the indices are the eigenvalues and the entries
    # are the Jordan structure corresponding to that eigenvalue.
    for pair in eigsAndSize do 
        if pair[1] in {indices(jordanStructure, 'nolist')} then 
            jordanStructure[pair[1]] := [op(jordanStructure[pair[1]]), pair[2]];
        else 
            jordanStructure[pair[1]] := [pair[2]];
        end if;
    end do;
    
    # Sort each Jordan structure list in decreasing order
    for item in indices(jordanStructure, 'nolist') do
        jordanStructure[item] := sort(jordanStructure[item], `>`);
    end do;
    
    return jordanStructure;
    
end proc;


# ----------------------------------------------------------------------- #
# eigsAndBlockSizeFromJCF                                                 #
#                                                                         #
# Get a list of pairs where each pair is [eig, blockSize]. The output     #
# list corresponds to the block structure of the input Jordan block       #
# matrix. The order in which the pairs appear in the output list matches  #
# that of the input matrix.                                               #
#                                                                         #
# INPUT                                                                   #
#   J ... Jordan block matrix                                             #
#                                                                         #
# OUTPUT                                                                  #
#   A list of pairs where each pair is [eig, blockSize] corresponding to  #
#   the Jordan blocks in the input matrix.                                #
# ----------------------------------------------------------------------- #
eigsAndBlockSizeFromJCF := proc(J::Matrix, $)::list([anything, posint]);

    local superDiagonal::list,
          blockLocation::list(truefalse),
             blockInter::list,
             splitBlock::list(list);

    # Get a list of the values on the first super diagonal
    superDiagonal:=convert(Diagonal(J,1),list); 

    blockLocation := map(proc (a) options operator, arrow; evalb(a = 1) end proc, superDiagonal);

    blockInter := Interleave(convert(Diagonal(J), list), blockLocation);

    splitBlock := [Split(`=`, blockInter, false)];

    splitBlock := map2(remove, proc (a) options operator, arrow; evalb(a = true) end proc, splitBlock);

    return map(proc (a) options operator, arrow; [a[1], nops(a)] end proc, splitBlock);
    
end proc;


# ----------------------------------------------------------------------- #
# jordanStructureToWeyrStructure                                          #
#                                                                         #
# Convert a Jordan structure to its Weyr structure.                       #
# Example:                                                                #
#   Jordan Structure = [3,2,2]                                            #
#   0 <- v1 <- v2 <- v3                                                   #
#   0 <- v5 <- v6                                                         #
#   0 <- v7 <- v8                                                         #
#                                                                         #
#   Transpose:                                                            #
#   0   0   0                                                             #
#   ^   ^   ^                                                             #
#   v1  v5  v7                                                            #
#   ^   ^   ^                                                             #
#   v2  v6  v8                                                            #
#   ^   ^   ^                                                             #
#   v4                                                                    #
#   => Weyr Structure = [3,3,1]                                           #
#                                                                         #
# INPUT                                                                   #
#   jordanStructure ... (list) The Jordan structure                       #
#                                                                         #
# OUTPUT                                                                  #
#   A list corresponding to the Weyr structure associated with the input  #
#   Jordan structure.                                                     #
#                                                                         #
# ASSUMPTIONS                                                             #
#   The entries in the input list are in decreasing order.                #
# ----------------------------------------------------------------------- #
jordanStructureToWeyrStructure := proc(jordanStructure::list(posint), $)::list(posint);
    
    local weyrStructure::list(posint), l::list(nonnegint); 
    
    # Copy jordanStructure so we can operate on it
    l := jordanStructure; 
    
    # Create an empty list for the Weyr structure
    weyrStructure := []; 
     
    while 0 < nops(l) do 
        weyrStructure := [op(weyrStructure), nops(l)]; 
        l := map(proc (a) options operator, arrow; a-1 end proc, l); 
        l := remove(proc (a) options operator, arrow; evalb(a = 0) end proc, l);
    end do; 
    
    return weyrStructure;
    
end proc;


# ----------------------------------------------------------------------- #
# weyrBlockMatrix                                                         #
#                                                                         #
# Build a Weyr block matrix with a single eigenvalue.                     #
#                                                                         #
# INPUT                                                                   #
#   eigVal .......... Eigenvalue                                          #
#   weyrStructure ... The Weyr structure for the output matrix            #
#                                                                         #
# OUTPUT                                                                  #
#   A block matrix where the diagonal blocks are eigVal*IdentityMatrix    #
#   of the size given in the WeyrStructure list. The super diagonal       #
#   blocks are identity matrices of the appropriate size.                 #
# ----------------------------------------------------------------------- #
weyrBlockMatrix := proc(eigVal, weyrStructure::(list(posint)), $)::Matrix;
    
    local superDiagBlockList::list(Matrix),
                           i::posint,
                        nRow::posint, 
                        nCol::posint, 
                           W::Matrix, 
                  matrixSize::posint,
              identityBlocks::Matrix,
               startBlockCol::posint,
                 endBlockCol::posint,
               startBlockRow::posint,
                 endBlockRow::posint;
    
    if nops(weyrStructure) = 1 then
        return eigVal*LinearAlgebra:-IdentityMatrix(weyrStructure[1]);
    end if;
    
    superDiagBlockList := [];
    
    for i to nops(weyrStructure)-1 do
        nRow := weyrStructure[i];
        nCol := weyrStructure[i+1];
        
        superDiagBlockList := [op(superDiagBlockList), LinearAlgebra:-IdentityMatrix(nRow, nCol)]
    end do;
    
    matrixSize := add(weyrStructure);
    W := Matrix(matrixSize);
    
    identityBlocks := LinearAlgebra:-DiagonalMatrix(superDiagBlockList);
    
    startBlockCol := weyrStructure[1]+1;
    endBlockCol := matrixSize;
    startBlockRow := 1;
    endBlockRow := matrixSize-weyrStructure[nops(weyrStructure)];
    
    W[startBlockRow .. endBlockRow, startBlockCol .. endBlockCol] := identityBlocks;
    
    for i to matrixSize do
        W[i, i] := eigVal;
    end do; 
    
    return W;

end proc;


# ----------------------------------------------------------------------- #
# sortJordanForm                                                          #
#                                                                         #
# Return the similarity transformation matrix Q such that J2 = Q^-1 J Q   #
# where both J and J2 are Jordan block matrices, but J2 has blocks        #
# corresponding to the same eigenvalues grouped together (i.e.            #
# sequentially along the diagonal) and the blocks for each eigenvalue are #
# in decreasing order by block size.                                      #
#                                                                         #
# INPUT                                                                   #
#   J ... Jordan block matrix                                             #
#                                                                         #
# OUTPUT                                                                  #
#   The similarity transformation matrix Q such that J2 = Q^-1 J Q where  #
#   J2 is a Jordan block matrix that is similar to the input matrix with  #
#   the Jordan blocks permuted along the diagonal such that block         #
#   corresponding to the same eigenvalues appear sequentially along the   #
#   diagonal. The grouped blocks for each eigenvalue appear in decreasing #
#   order of block size along the diagonal.                               #
# ----------------------------------------------------------------------- #
sortJordanForm := proc(J::Matrix, $)
    
    local n::posint,
          Q::Matrix,
          Q2::Matrix,
    eigVals::list,
        eig,
          i::posint,
          permutation::list(posint),
          J2::Matrix,
          blockList, startIndex, endIndex, permutationList, J3;
    
    n := RowDimension(J);
    
    if n = 1 then
        return Matrix([1]);
    end if;
    
    # A set of the unique eigenvalues in the matrix
    eigVals := MakeUnique(convert(Diagonal(J), list));
    
    # 1. Group blocks corresponding to the same eigenvalue
    permutation := [];
    for eig in eigVals do
        for i to n do
            if J[i,i] = eig then
                permutation := [op(permutation), i];
            end if;
        end do;
    end do;
    
    Q := permutationMatrix(permutation);
    Q := MatrixInverse(Q);
    
    # 2. Order blocks in decreasing order
    J2 := MatrixInverse(Q).J.Q;
    
    ASSERT(IsSimilar(J, J2), "J2 must be similar to J");
    
    # Get the eigenvalues (in the same order they appear in J2)
    eigVals := MakeUnique(convert(Diagonal(J2), list));
    
    # Split into Jordan block matrices for each eigenvalue
    blockList := [];
    
    startIndex := 1;
    endIndex := 1;
    for eig in eigVals do
        while endIndex <= n and J2[endIndex, endIndex] = eig do
            endIndex := endIndex + 1;
        end do;
        
        blockList := [op(blockList), J2[startIndex..(endIndex-1),startIndex..(endIndex-1)]];
        
        startIndex := endIndex;
    end do;
    
    # Sort each Jordan block matrix in blockList
    permutationList := [];
    for J3 in blockList do
        eig := J3[1,1];
        permutationList := [op(permutationList), sortJordanBlock(J3)];
    end do;
    
    Q2 := DiagonalMatrix(permutationList);
    
    return Q.Q2;
    
end proc;


# ----------------------------------------------------------------------- #
# sortJordanBlock                                                         #
#                                                                         #
# Return the permutation matrix Q such that J2 = Q^-1 J Q where J is a    #
# Jordan block matrix with one eigenvalue and J2 is also a Jordan block   #
# matrix where the blocks appear in decreasing order along the diagonal.  #
#                                                                         #
# INPUT                                                                   #
#   J ... Jordan block matrix corresponding to a single eigenvalue        #
#                                                                         #
# OUTPUT                                                                  #
#   The similarity transformation matrix Q such that J2 = Q^-1 J Q where  #
#   J2 is a Jordan block matrix that is similar to the input matrix with  #
#   the Jordan blocks permuted along the diagonal such that the blocks    #
#   are ordered by decreasing size.                                       #
# ----------------------------------------------------------------------- #
sortJordanBlock := proc(J::Matrix, $)::Matrix;
    
    local blockStructure::list([anything, posint]),
              blockSizes::list(posint), 
                       i::posint, 
       originalPositions::list(nonnegint), 
                       P::list(posint), 
            newPositions::list(nonnegint), 
         currentPosition::posint, 
             endRowIndex::list(nonnegint), 
             endColIndex::list(nonnegint), 
                       Q::Matrix;
    
    if RowDimension(J) = 1 then
        return Matrix([1]);
    end if;
    
    blockStructure := eigsAndBlockSizeFromJCF(J);
    
    blockSizes := map(proc (a) options operator, arrow; a[2] end proc, blockStructure);
    
    originalPositions := [0$nops(blockSizes)]; 
    originalPositions[1] := 1;
    for i from 2 to nops(blockSizes) do 
        originalPositions[i] := originalPositions[i-1]+blockSizes[i-1];
    end do;
    
    P := sort(blockSizes, proc (x, y) options operator, arrow; y < x end proc, 'output' = 'permutation');
    
    newPositions := [0$nops(blockSizes)]; 
    currentPosition := 1; 
    for i to nops(blockSizes) do 
        newPositions[P[i]] := currentPosition; 
        currentPosition := currentPosition+blockSizes[P[i]];
    end do; 
    
    endRowIndex := zip(`-`, zip(`+`, blockSizes, newPositions), 1);
    endColIndex := zip(`-`, zip(`+`, blockSizes, originalPositions), 1);
    
    Q := Matrix(add(blockSizes)); 
    for i to nops(blockSizes) do 
        Q[newPositions[i] .. endRowIndex[i], originalPositions[i] .. endColIndex[i]] := IdentityMatrix(blockSizes[i], 'compact' = false);
    end do;
    
    Q := MatrixInverse(Q);
    
    return Q;
    
end proc;


# ----------------------------------------------------------------------- #
# JCF_to_WCF_Transformation_One_Eig                                       #
#                                                                         #
# Compute the transformation matrix from a matrix in JCF to WCF where the #
# input matrix has a single eigenvalue.                                   #
#                                                                         #
# INPUT                                                                   #
#   J ... Jordan block matrix with one eigenvalue.                        #
#                                                                         #
# OUTPUT                                                                  #
#   A permutation matrix as described above.                              #
# ----------------------------------------------------------------------- #
JCF_to_WCF_Transformation_One_Eig := proc(J::Matrix, $)
    
    local n, jordanStructure, A, count, i, j;
    
    n := RowDimension(J);
    
    if n = 1 then
        return Matrix([1]);
    end if;
    
    # Get the Weyr structure
    jordanStructure := op(entries(getJordanStructure(J)));
    
    # Make a matrix with the Jordan structure
    A := Matrix(nops(jordanStructure), max(jordanStructure));
    count := 1;
    for i to nops(jordanStructure) do
        for j to jordanStructure[i] do
            A[i,j] := count;
            count := count + 1;
        end do;
    end do;
    
    A := convert(A, list);
    A := remove(proc (a) options operator, arrow; a = 0 end proc, A);
    
    return permutationMatrix(A);
    
end proc;


# ----------------------------------------------------------------------- #
# permutationMatrix                                                       #
#                                                                         #
# Return a square matrix with dimension equal to the number of elements   #
# in the input list. The entries of the matrix are all zeros except for a #
# 1 in the ith column on the jth row where j is the index in the input    #
# list and i is the jth value.                                            #
#                                                                         #
# INPUT                                                                   #
#   P ... A list of positive integers where no value appears more than    #
#         once and it contains all integers values in {1..max(P)}         #
#                                                                         #
# OUTPUT                                                                  #
#   A permutation matrix as described above.                              #
# ----------------------------------------------------------------------- #
permutationMatrix := proc(P::list(posint), $)::Matrix;
    
    local A::Matrix, i::posint;
    
    A := Matrix(nops(P));
    
    for i to nops(P) do
        A[i, P[i]] := 1;
    end do;
    
    return A;
    
end proc;

end module: