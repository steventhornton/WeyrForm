# ======================================================================= #
# ======================================================================= #
#                                                                         #
# test.mpl                                                                #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 13/2016                                                #
#                                                                         #
# Run test for the WeyrForm function and all sub-procedures.              #
#                                                                         #
# CALLING SEQUENCE                                                        #
#   maple -q test.mpl                                                     #
#                                                                         #
# LICENSE                                                                 #
#   This program is free software: you can redistribute it and/or modify  #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation, either version 3 of the License, or     #
#   any later version.                                                    #
#                                                                         #
#   This program is distributed in the hope that it will be useful,       #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
#   You should have received a copy of the GNU General Public License     #
#   along with this program.  If not, see http://www.gnu.org/licenses/.   #
# ======================================================================= #
# ======================================================================= #
read("WeyrForm.mpl"):
kernelopts('opaquemodules' = false):
kernelopts('assertlevel' = 2):
infolevel['WeyrForm'] := 2:
with(LinearAlgebra):

# ----------------------------------------------------------------------- #
# WeyrForm                                                                #
# ----------------------------------------------------------------------- #
test101 := proc()
    
    local A, QQ, WW, W_correct, Q_correct;
    
    A := Matrix([[-1, 0, 0, 2, 1, 0, 0], 
                 [-1, 0, 0, 1, 1, 0, 0], 
                 [6, -2, -2, 4, 2, -2, 4], 
                 [-2, 1, 1, 0, 0, 1, -1], 
                 [3, -1, -1, 2, 1, -1, 2], 
                 [-5, 2, 2, -5, -3, 2, -4], 
                 [2, 0, 0, -4, -2, 0, 0]]);
    W_correct := Matrix([[0, 0, 0, 1, 0, 0, 0], 
                         [0, 0, 0, 0, 1, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 1, 0], 
                         [0, 0, 0, 0, 0, 0, 1], 
                         [0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0]]);
    Q_correct := Matrix([[1, 1, 1/4, 5/12, -7/12, 1/48, 49/48], 
                         [-1/3, 5/3, 1/4, 5/12, -7/12, 3/16, 3/16], 
                         [2, 2, -3/2, -5/2, 7/2, -9/8, -9/8], 
                         [0, 0, 0, 4/3, -2/3, 0, 0], 
                         [1, 1, 1/4, -5/4, 7/4, 7/16, 7/16], 
                         [-5/3, -11/3, 5/4, 25/12, -35/12, 31/16, 31/16], 
                         [-2, -2, -1/2, -5/6, 7/6, -3/8, -3/8]]);
    
    QQ, WW := WeyrForm(A, 'output'=['Q', 'W']);
    
    if not Equal(WW, W_correct, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    if not Equal(QQ, Q_correct, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
    
end proc:


test102 := proc()

    local A, WW, QQ;
    
    to 10 do
        A := RandomMatrix(3);
        
        WW, QQ := WeyrForm(A, 'output'=['W', 'Q']);
        
        if not Equal(Matrix(3), map(simplify, WW - MatrixInverse(QQ).A.QQ)) then
            printf("\n");
            error "Wrong result";
        end if;
        
    end do;
    
    printf("Passed!\n");
    
end proc:


test103 := proc()

    local A, QQ, WW, W_correct, Q_correct;

    A := Matrix([[0]]);
    W_correct := Matrix([[0]]);
    Q_correct := Matrix([[1]]);

    QQ, WW := WeyrForm(A, 'output'=['Q', 'W']);

    if not Equal(WW, W_correct, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    if not Equal(QQ, Q_correct, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
    
end proc:


test104 := proc()
    
    local J, WW, QQ, Z;
    
    J := Matrix([[5, 0, 0, 0, 0, 0, 0], 
                 [0, 2, 1, 0, 0, 0, 0], 
                 [0, 0, 2, 0, 0, 0, 0], 
                 [0, 0, 0, 5, 1, 0, 0], 
                 [0, 0, 0, 0, 5, 0, 0], 
                 [0, 0, 0, 0, 0, 2, 0], 
                 [0, 0, 0, 0, 0, 0, 5]]);
    
    WW, QQ := WeyrForm(J, 'output'=['W', 'Q']);
    
    Z :=  MatrixInverse(QQ).J.QQ - WW;
    
    Z := map(simplify, Z);
    
    if not Equal(Z, Matrix(7), 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
end proc:


test105 := proc()

    local A, WW, QQ, Z;

    A := Matrix([[2, 1, 0, 0, 0, 0, 0, 0],
                 [0, 2, 1, 0, 0, 0, 0, 0],
                 [0, 0, 2, 0, 0, 0, 0, 0],
                 [0, 0, 0, 2, 1, 0, 0, 0],
                 [0, 0, 0, 0, 2, 0, 0, 0],
                 [0, 0, 0, 0, 0, 2, 1, 0],
                 [0, 0, 0, 0, 0, 0, 2, 0],
                 [0, 0, 0, 0, 0, 0, 0, 2]]);

    WW, QQ := WeyrForm(A, 'output'=['W', 'Q']);

    Z :=  MatrixInverse(QQ).A.QQ - WW;

    Z := map(simplify, Z);

    if not Equal(Z, Matrix(8), 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
end proc:


test106 := proc()

    local A, WW, QQ, Z;

    A := Matrix([[0, -3, 1, 2],
                 [-2, 1, -1, 2],
                 [-2, 1, -1, 2],
                 [-2, -3, 1, 4]]);

    WW, QQ := WeyrForm(A, 'output'=['W', 'Q']);

    Z :=  MatrixInverse(QQ).A.QQ - WW;

    Z := map(simplify, Z);

    if not Equal(Z, Matrix(4), 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
end proc:


# ----------------------------------------------------------------------- #
# WeyrForm:-implementation_W                                              #
# ----------------------------------------------------------------------- #
#test201 := proc()

# ----------------------------------------------------------------------- #
# WeyrForm:-implementation_WQ                                             #
# ----------------------------------------------------------------------- #
#test301 := proc()

# ----------------------------------------------------------------------- #
# WeyrForm:-getJordanStructure                                            #
# ----------------------------------------------------------------------- #
test401 := proc()
    
    local J, jStruct, jStructCorrect;
    
    J := Matrix(6); J[1,2] := 1;
    
    jStructCorrect := table([0=[2, 1, 1, 1, 1]]);
    
    jStruct := WeyrForm:-getJordanStructure(J);
    
    if not verify(jStructCorrect, jStruct, 'table') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
    
end proc:


test402 := proc()

    local J, jStruct, jStructCorrect;

    J := JordanBlockMatrix([[1,2], [2,2], [1,2], [1,1], [2,1]]);

    jStructCorrect := table([1=[2, 2, 1], 2=[2, 1]]);

    jStruct := WeyrForm:-getJordanStructure(J);

    if not verify(jStructCorrect, jStruct, 'table') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test403 := proc()

    local J, jStruct, jStructCorrect;

    J := DiagonalMatrix([0, 1, 2, 3, 4]);

    jStructCorrect := table([0=[1], 1=[1], 2=[1], 3=[1], 4=[1]]);

    jStruct := WeyrForm:-getJordanStructure(J);

    if not verify(jStructCorrect, jStruct, 'table') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test404 := proc()

    local J, jStruct, jStructCorrect;

    J := JordanBlockMatrix([[-1, 18]]);

    jStructCorrect := table([-1=[18]]);

    jStruct := WeyrForm:-getJordanStructure(J);

    if not verify(jStructCorrect, jStruct, 'table') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


# ----------------------------------------------------------------------- #
# WeyrForm:-eigsAndBlockSizeFromJCF                                       #
# ----------------------------------------------------------------------- #
test501 := proc()
    
    local J, eigsAndSize, correct;
    
    correct := [[0, 2], [0, 1], [0, 1], [0, 1], [0, 1]];
    
    J := JordanBlockMatrix(correct);
    
    eigsAndSize := WeyrForm:-eigsAndBlockSizeFromJCF(J);
    
    if not evalb(correct = eigsAndSize) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
    
end proc:


test502 := proc()

    local J, eigsAndSize, correct;

    correct := [[1, 2], [2, 2], [1, 2], [1, 1], [2, 1]];

    J := JordanBlockMatrix(correct);

    eigsAndSize := WeyrForm:-eigsAndBlockSizeFromJCF(J);

    if not evalb(correct = eigsAndSize) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test503 := proc()

    local J, eigsAndSize, correct;

    correct := [[-1, 1], [0, 1], [1, 1], [2, 1], [3, 1], [4, 1]];

    J := JordanBlockMatrix(correct);

    eigsAndSize := WeyrForm:-eigsAndBlockSizeFromJCF(J);

    if not evalb(correct = eigsAndSize) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test504 := proc()

    local J, eigsAndSize, correct;

    correct := [[-1, 18]];

    J := JordanBlockMatrix(correct);

    eigsAndSize := WeyrForm:-eigsAndBlockSizeFromJCF(J);

    if not evalb(correct = eigsAndSize) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test505 := proc()

    local J, eigsAndSize, correct;

    correct := [[-1, 1]];

    J := JordanBlockMatrix(correct);

    eigsAndSize := WeyrForm:-eigsAndBlockSizeFromJCF(J);

    if not evalb(correct = eigsAndSize) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


# ----------------------------------------------------------------------- #
# WeyrForm:-jordanStructureToWeyrStructure                                #
# ----------------------------------------------------------------------- #
test601 := proc()
    
    local jStruct, wStruct, wStructCorrect;
    
    jStruct := [3, 2, 2];
    wStructCorrect := [3, 3, 1];
    wStruct := WeyrForm:-jordanStructureToWeyrStructure(jStruct);
    
    if not evalb(wStruct = wStructCorrect) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test602 := proc()

    local jStruct, wStruct, wStructCorrect;

    jStruct := [1];
    wStructCorrect := [1];
    wStruct := WeyrForm:-jordanStructureToWeyrStructure(jStruct);

    if not evalb(wStruct = wStructCorrect) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test603 := proc()

    local jStruct, wStruct, wStructCorrect;

    jStruct := [1, 1, 1, 1, 1];
    wStructCorrect := [5];
    wStruct := WeyrForm:-jordanStructureToWeyrStructure(jStruct);

    if not evalb(wStruct = wStructCorrect) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test604 := proc()

    local jStruct, wStruct, wStructCorrect;

    jStruct := [7];
    wStructCorrect := [1, 1, 1, 1, 1, 1, 1];
    wStruct := WeyrForm:-jordanStructureToWeyrStructure(jStruct);

    if not evalb(wStruct = wStructCorrect) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


# ----------------------------------------------------------------------- #
# WeyrForm:-weyrBlockMatrix                                               #
# ----------------------------------------------------------------------- #
test701 := proc()
    
    local W_correct, W, eigVal, weyrStructure;
    
    eigVal := 0;
    weyrStructure := [1];
    W_correct := Matrix(1, 1, [[0]]);
    
    W := WeyrForm:-weyrBlockMatrix(eigVal, weyrStructure);
    
    if not Equal(W, W_correct, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
    
end proc:


test702 := proc()

    local W_correct, W, eigVal, weyrStructure;
    
    eigVal := 3;
    weyrStructure := [1, 1, 1, 1, 1];
    W_correct := Matrix(5, 5, [[3, 1, 0, 0, 0], 
                               [0, 3, 1, 0, 0], 
                               [0, 0, 3, 1, 0], 
                               [0, 0, 0, 3, 1], 
                               [0, 0, 0, 0, 3]]);
    W := WeyrForm:-weyrBlockMatrix(eigVal, weyrStructure);
    
    if not Equal(W, W_correct, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
    
end proc:


test703 := proc()

    local W_correct, W, eigVal, weyrStructure;
    
    eigVal := -1;
    weyrStructure := [4, 3, 2, 1];
    W_correct := Matrix(10, 10, [[-1, 0, 0, 0, 1, 0, 0, 0, 0, 0], 
                                 [0, -1, 0, 0, 0, 1, 0, 0, 0, 0], 
                                 [0, 0, -1, 0, 0, 0, 1, 0, 0, 0], 
                                 [0, 0, 0, -1, 0, 0, 0, 0, 0, 0], 
                                 [0, 0, 0, 0, -1, 0, 0, 1, 0, 0], 
                                 [0, 0, 0, 0, 0, -1, 0, 0, 1, 0], 
                                 [0, 0, 0, 0, 0, 0, -1, 0, 0, 0], 
                                 [0, 0, 0, 0, 0, 0, 0, -1, 0, 1], 
                                 [0, 0, 0, 0, 0, 0, 0, 0, -1, 0], 
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, -1]]);
    
    W := WeyrForm:-weyrBlockMatrix(eigVal, weyrStructure);
    
    if not Equal(W, W_correct, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test704 := proc()

    local W_correct, W, eigVal, weyrStructure;
    
    eigVal := 3;
    weyrStructure := [4, 4, 3, 3, 2, 2, 1];
    W_correct := Matrix([[3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 1, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 1], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0], 
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3]]);
    
    W := WeyrForm:-weyrBlockMatrix(eigVal, weyrStructure);
    
    if not Equal(W, W_correct, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
    
end proc:


test705 := proc()

    local W_correct, W, eigVal, weyrStructure;
    
    eigVal := 8;
    weyrStructure := [4];
    W_correct := Matrix(4, 4, [[8, 0, 0, 0], 
                               [0, 8, 0, 0], 
                               [0, 0, 8, 0], 
                               [0, 0, 0, 8]]);
    
    W := WeyrForm:-weyrBlockMatrix(eigVal, weyrStructure);
    
    if not Equal(W, W_correct, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
    
end proc:


# ----------------------------------------------------------------------- #
# WeyrForm:-sortJordanForm                                                #
# ----------------------------------------------------------------------- #
test801 := proc()
    
    local J, J_Correct, Q, Z;
    
    J := JordanBlockMatrix([[0, 1],[0, 3],[0, 2],[0, 1], [0, 1]]);
    J_Correct := JordanBlockMatrix([[0, 3], [0, 2], [0, 1], [0, 1], [0, 1]]);
    
    Q := WeyrForm:-sortJordanForm(J);
    
    Z := MatrixInverse(Q).J.Q;
    Z := map(simplify, Z);
    
    if not Equal(Z, J_Correct, 'compare'='entries') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
    
end proc:


test802 := proc()

    local J, J_Correct, Q, Z;

    J := JordanBlockMatrix([[1, 2], [2, 2], [1, 2], [1, 1], [2, 1]]);
    J_Correct := JordanBlockMatrix([[1, 2], [1, 2], [1, 1], [2, 2], [2, 1]]);

    Q := WeyrForm:-sortJordanForm(J);

    Z := MatrixInverse(Q).J.Q;
    Z := map(simplify, Z);

    if not IsSimilar(Z, J_Correct) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test803 := proc()

    local J, Q, Z;

    J := JordanBlockMatrix([[-1, 1], [0, 1], [1, 1], [2, 1], [3, 1]]);

    Q := WeyrForm:-sortJordanForm(J);

    Z := MatrixInverse(Q).J.Q;
    Z := map(simplify, Z);

    if not IsSimilar(Z, J) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test804 := proc()

    local J, Q, Z;

    J := JordanBlockMatrix([[-1, 1]]);

    Q := WeyrForm:-sortJordanForm(J);

    Z := MatrixInverse(Q).J.Q;
    Z := map(simplify, Z);

    if not IsSimilar(Z, J) then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


# ----------------------------------------------------------------------- #
# WeyrForm:-sortJordanBlock                                               #
# ----------------------------------------------------------------------- #
test901 := proc()

    local J, J_Correct, Q, Z;

    J := JordanBlockMatrix([[0, 1],[0, 3],[0, 2],[0, 1], [0, 1]]);
    J_Correct := JordanBlockMatrix([[0, 3], [0, 2], [0, 1], [0, 1], [0, 1]]);

    Q := WeyrForm:-sortJordanForm(J);

    Z := MatrixInverse(Q).J.Q;
    Z := map(simplify, Z);

    if not Equal(Z, J_Correct, 'compare'='entries') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


test902 := proc()

    local J, Q, Z;

    J := JordanBlockMatrix([[0, 1]]);

    Q := WeyrForm:-sortJordanForm(J);

    Z := MatrixInverse(Q).J.Q;
    Z := map(simplify, Z);

    if not Equal(Z, J, 'compare'='entries') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


# ----------------------------------------------------------------------- #
# WeyrForm:-JCF_to_WCF_Transformation_One_Eig                             #
# ----------------------------------------------------------------------- #
test1001 := proc()
    
    local J, W, Q, Z;
    
    J := Matrix([[10]]);
    W := Matrix([[10]]);
    
    Q := WeyrForm:-JCF_to_WCF_Transformation_One_Eig(J);
    
    Z := Q.J.MatrixInverse(Q);
    
    Z := map(simplify, Z);
    
    if not Equal(Z, W, 'compare'='entries') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
    
end proc:


test1002 := proc()

    local J, W, Q, Z;

    J := JordanBlockMatrix([[0, 2], [0, 1], [0, 1], [0, 1], [0, 1]]);
    W := Matrix(6); W[1,6] := 1;

    Q := WeyrForm:-JCF_to_WCF_Transformation_One_Eig(J);

    Z := Q.J.MatrixInverse(Q);

    Z := map(simplify, Z);
    
    if not Equal(Z, W, 'compare'='entries') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");

end proc:


# ----------------------------------------------------------------------- #
# WeyrForm:-permutationMatrix                                             #
# ----------------------------------------------------------------------- #
test1101 := proc()
    local A, P, l;
    l := [1, 2, 3];
    A := Matrix([[1,0,0],[0,1,0],[0,0,1]]);
    P := WeyrForm:-permutationMatrix(l);
    if not Equal(A, P, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
end proc:


test1102 := proc()
    local A, P, l;
    l := [1, 5, 6, 12, 7, 2, 3, 4, 10, 9, 11, 8];
    A := Matrix([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                 [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                 [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                 [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]]);
    P := WeyrForm:-permutationMatrix(l);
    if not Equal(A, P, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
end proc:


test1103 := proc()
    local A, P, l;
    l := [1];
    A := Matrix([[1]]);
    P := WeyrForm:-permutationMatrix(l);
    if not Equal(A, P, 'compare'='all') then
        printf("\n");
        error "Wrong result";
    end if;
    printf("Passed!\n");
end proc:


# ----------------------------------------------------------------------- #
# Run tests                                                               #
# ----------------------------------------------------------------------- #

printf("WeyrForm\n");

printf("\tTest 101: ");
test101();

printf("\tTest 102: ");
test102();

printf("\tTest 103: ");
test103();

printf("\tTest 104: ");
test104();

printf("\tTest 105: ");
test105();

printf("\tTest 106: ");
test106();


printf("WeyrForm:-getJordanStructure\n");

printf("\tTest 401: ");
test401();

printf("\tTest 402: ");
test402();

printf("\tTest 403: ");
test403();

printf("\tTest 404: ");
test404();


printf("WeyrForm:-eigsAndBlockSizeFromJCF  \n");

printf("\tTest 501: ");
test501();

printf("\tTest 502: ");
test502();

printf("\tTest 503: ");
test503();

printf("\tTest 504: ");
test504();

printf("\tTest 505: ");
test505();


printf("WeyrForm:-jordanStructureToWeyrStructure\n");

printf("\tTest 601: ");
test601();

printf("\tTest 602: ");
test602();

printf("\tTest 603: ");
test603();

printf("\tTest 604: ");
test604();


printf("WeyrForm:-weyrBlockMatrix\n");

printf("\tTest 701: ");
test701();

printf("\tTest 702: ");
test702();

printf("\tTest 703: ");
test703();

printf("\tTest 704: ");
test704();

printf("\tTest 705: ");
test705();


printf("WeyrForm:-sortJordanForm\n");

printf("\tTest 801: ");
test801();

printf("\tTest 802: ");
test802();

printf("\tTest 803: ");
test803();

printf("\tTest 804: ");
test804();


printf("WeyrForm:-sortJordanBlock\n");

printf("\tTest 901: ");
test901();

printf("\tTest 902: ");
test902();


printf("WeyrForm:-JCF_to_WCF_Transformation_One_Eig \n");

printf("\tTest 1001: ");
test1001();

printf("\tTest 1002: ");
test1002();


printf("WeyrForm:-permutationMatrix\n");

printf("\tTest 1101: ");
test1101();

printf("\tTest 1102: ");
test1102();

printf("\tTest 1103: ");
test1103();