# ======================================================================= #
# ======================================================================= #
#                                                                         #
# test.mpl                                                                #
#                                                                         #
# AUTHOR .... Steven E. Thornton                                          #
# EMAIL ..... sthornt7@uwo.ca                                             #
# UPDATED ... Oct. 8/2016                                                 #
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

    local A, WW, QQ, i;
    
    for i to 10 do
        A := RandomMatrix(3);
<<<<<<< HEAD
=======
        printf("\t%d of 100\n", i);
>>>>>>> 325c59848b0a51c8f64969ff97e55a7c1734f811
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
#test401 := proc()

# ----------------------------------------------------------------------- #
# WeyrForm:-eigsAndBlockSizeFromJCF                                       #
# ----------------------------------------------------------------------- #
#test501 := proc()

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
    W_correct := Matrix(19, 19, [[3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3]]);
    
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
#test801 := proc()

# ----------------------------------------------------------------------- #
# WeyrForm:-sortJordanBlock                                               #
# ----------------------------------------------------------------------- #
#test901 := proc()

# ----------------------------------------------------------------------- #
# WeyrForm:-JCF_to_WCF_Transformation_One_Eig                             #
# ----------------------------------------------------------------------- #
#test1001 := proc()

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


printf("WeyrForm:-permutationMatrix\n");

printf("\tTest 1101: ");
test1101();

printf("\tTest 1102: ");
test1102();

printf("\tTest 1103: ");
test1103();