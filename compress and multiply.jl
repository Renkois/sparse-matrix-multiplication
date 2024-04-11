using SparseArrays
function matrix_multiply(A::SparseMatrixCSC, B::SparseMatrixCSC)    #用CSC的方式进行优化
	sidelength = size(A,1)
	Anon_zero = A.nzval; Aposition = A.rowval; Afirst = A.colptr    #分别获取A中的非零数、每个非零所在数的行数、每一列第一个非零数的位置
    Bnon_zero = B.nzval; Bposition = B.rowval; Bfirst = B.colptr
	result = zeros(sidelength,sidelength)
    @inbounds for i = 1:sidelength
        if Afirst[i+1] - Afirst[i] != 0
            @inbounds for j = Afirst[i]:Afirst[i+1] - 1
                if Bfirst[Aposition[j]+1] - Bfirst[Aposition[j]] != 0
                    @inbounds for k = Bfirst[Aposition[j]]:Bfirst[Aposition[j]+1] - 1
                        result[Bposition[k],i] += Anon_zero[j] * Bnon_zero[k]
                    end
                end
            end
        end
    end
    return result
end
