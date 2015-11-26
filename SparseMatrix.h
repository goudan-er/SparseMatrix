#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
using namespace std;

namespace AHaTinyMath
{
	template <class T>
	class SparseMatrix {
	public:
		typedef vector< unordered_map<size_t, size_t> > map_mat_t;
		typedef unordered_map<size_t, size_t> map_col_t;
		typedef typename map_col_t::iterator map_col_iter;

		typedef vector< map<size_t, T> > mat_t;
		typedef map<size_t, T> col_t;
		typedef typename col_t::iterator col_iter;

		typedef vector< pair<size_t, size_t> > outline_t;
		//!+ Notice, tuple is the new feature of C++ 11
		typedef vector < tuple<size_t, size_t, T> > outline_with_value_t;
		
		// constructor
		SparseMatrix(size_t _numRows, size_t _numColumns)
			: numRows(_numRows), numColumns(_numColumns)
		{
			Allocate();
		}
		SparseMatrix(const SparseMatrix & source)
		{
			numRows = source.GetNumRows();
			numColumns = source.GetNumColumns();
			map_mat.resize(numRows);

			rowLength.resize(numRows);
			columnIndices.resize(numRows);
			columnEntries.resize(numRows);
			for (size_t i = 0, j; i < numRows; ++i) { // i: row, j: column
				rowLength[i] = source.rowLength[i];
				columnIndices[i].resize(rowLength[i]);
				columnEntries[i].resize(rowLength[i]);
				for (size_t k = 0; k < rowLength[i]; ++k) {
					columnIndices[i][k] = source.columnIndices[i][k];
					columnEntries[i][k] = source.columnEntries[i][k];
					j = columnIndices[i][k];
					map_mat[i][j] = k;
				}
			}

		}
		// destructor
		~SparseMatrix()
		{
			for (size_t i = 0; i < numRows; ++i) {
				map_mat[i].clear();
			}
			map_mat.clear();

			for (size_t i = 0; i < numRows; ++i) {
				columnIndices[i].clear();
				columnEntries[i].clear();
			}
			columnIndices.clear();
			columnEntries.clear();

			rowLength.clear();
		}

		size_t GetNumRows() const { return numRows; }
		size_t GetNumColumns() const { return numColumns; }

		size_t* GetRowLength() const { return &rowLength[0]; }
		vector<size_t>* GetColumnIndices() const { return &columnIndices[0]; }
		vector<size_t>* GetColumnEntries() const { return &columnEntries[0]; }

		// give all non-zero entries's position (i, j) of the SparseMatrix
		// and all non-zero entries initialized by value
		void InitFromOutline(outline_t & outline, const T & value = 0)
		{
			mat_t mat;
			// allocate mat
			mat.clear();
			col_t emptyRow;
			for (size_t i = 0; i < numRows; ++i) {
				mat.push_back(emptyRow);
			}
			// first, based outline & value, generate mat and it storaged by a map
			outline_t::iterator	pBegin = outline.begin(), pEnd = outline.end();
			for (outline_t::iterator it = pBegin; it < pEnd; ++it) {
				size_t i = it->first, j = it->second;
				mat[i][j] = value;
			}
			// second, to accelerate access mat[i][j], use vector storage column entries and indices
			// use map_mat to map the column to index
			for (size_t i = 0, j, k; i < numRows; ++i) {
				rowLength[i] = mat[i].size();
				k = 0; // the index of entry in each row
				for (col_iter it = mat[i].begin(); it != mat[i].end(); ++it) {
					columnIndices[i].push_back(it->first); // j
					columnEntries[i].push_back(it->second); // entry
					map_mat[i][it->first] = k++;
				}
				// deallocate column entries
				mat[i].clear();
			}
			mat.clear();
		}
		
		// the outline is tuple<i, j, value>
		// tuple is the new feature of C++ 11, you should guarantee your compilers support it
		void InitFromOutline(outline_with_value_t & outline)
		{
			mat_t mat;
			// allocate mat
			mat.clear();
			col_t emptyRow;
			for (size_t i = 0; i < numRows; ++i) {
				mat.push_back(emptyRow);
			}
			// first, based outline & value, generate mat and it storaged by a map
			outline_with_value_t::iterator	pBegin = outline.begin(), pEnd = outline.end();
			for (outline_with_value_t::iterator it = pBegin; it < pEnd; ++it) {
				size_t i = get<0>(*it), j = get<1>(*it);
				mat[i][j] = get<2>(*it);
			}
			// second, to accelerate access mat[i][j], use vector storage column entries and indices
			// use map_mat to map the column to index
			for (size_t i = 0, j, k; i < numRows; ++i) {
				rowLength[i] = mat[i].size();
				k = 0; // the index of entry in each row
				for (col_iter it = mat[i].begin(); it != mat[i].end(); ++it) {
					columnIndices[i].push_back(it->first); // j
					columnEntries[i].push_back(it->second); // entry
					map_mat[i][it->first] = k++;
				}
				// deallocate column entries
				mat[i].clear();
			}
			mat.clear();
		}
		
		// Intialized to diagonal matrix, the diagonal entries all are value
		// By default, value = 1,  i.e. the matrix is identity matrix
		void InitToDiagonal(const T & value = 1)
		{
			outline_t diagOutline;
			for (size_t i = 0; i < numRows; ++i) {
				diagOutline.push_back(make_pair(i, i));
			}
			InitFromOutline(diagOutline, value);
			diagOutline.clear();
		}
		
		// Intialized from an array data
		void InitFromDataPtr(T * data)
		{
			for (size_t i = 0; i < numRows; ++i) {
				for (size_t j = 0; j < numColumns; j++) {
					if (data[i*numColumns + j] != 0) {
						columnIndices[i].push_back(j);
						columnEntries[i].push_back(data[i*numColumns + j]);
						map_mat[i][j] = rowLength[i]++;
					}
				}
			}
		}

		//TODO 读入矩阵文件初始化
		void InitFromFile(char* filename);

		// Set all entries to zero
		void ResetZero()
		{
			for (size_t i = 0; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					columnEntries[i][k] = 0;
				}
			}
		}

		// 检查是否有非零元(i, j)
		// true: 存在, 并返回j在第i行的存储序号k
		// false: 不存在
		bool CheckEntry(size_t i, size_t j, size_t & k)
		{
			if (i >= numRows || map_mat[i].find(j) == map_mat[i].end()) {
				//std::cerr << "Error, the Sparse Matrix don't have (" << i << ',' << j << ") entry." << std::endl;
				return false;
			}
			k = map_mat[i][j];
			return true;
		}

		// Check matrix whether it is sysmetric
		bool CheckSysmetric()
		{
			if (numRows != numColumns) return false;
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					if (CheckEntry(j, i) && mat[i][j] == mat[j][i]) continue;
					else return false;
				}
			}
			return true;
		}

		// Set/Get/Add Entry, The Sparse Matrix must have the noo-zero entry (i, j)
		void SetEntry(size_t i, size_t j, const T & value = 0)
		{
			size_t k = 0;
			CheckEntry(i, j, k);
			//assert(CheckEntry(i, j, k));
			columnEntries[i][k] = value;
		}
		T GetEntry(size_t i, size_t j)
		{
			size_t k = 0;
			CheckEntry(i, j, k);
			//assert(CheckEntry(i, j, k));
			return columnEntries[i][k];
		}
		void AddEntry(size_t i, size_t j, const T & value)
		{
			size_t k = 0;
			CheckEntry(i, j, k);
			//assert(CheckEntry(i, j, k));
			columnEntries[i][k] += value;
		}
		
		//! 添加新的非零元，需要更新 rowLength & columnIndices & columnEntries，相对耗时
		void AddNewEntry(size_t i, size_t j, const T & value)
		{
			size_t k = 0;
			if (CheckEntry(i, j, k)) {
				columnEntries[i][k] += value;
			}
			else {
				rowLength[i] += 1;
				// 为了维护 columnIndices 的从小到大顺序, 搜索到对应的位置 *it
				vector<size_t>::iterator it = lower_bound(columnIndices[i].begin(), columnIndices[i].end(), j);
				// 插入 column j
				vector<size_t>::iterator pos = columnIndices[i].insert(it, j);
				// 插入 j entry, 同时建立 j-index 的索引
				vector<T>::iterator entryIt = columnEntries[i].begin();
				entryIt += *pos - 1;
				columnEntries[i].insert(entryIt, value);
				map_mat[i][j] = *pos - 1;
			}
		}

		//+ Matrix Algebra, 重载运算符
		//+ Notice: 运算的矩阵与mat有一样的pattern，即非零元相对应
		SparseMatrix operator+(const SparseMatrix & mat2) const
		{
			SparseMatrix result(*this);
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					result.columnEntries[i][k] += mat2.columnEntries[i][k];
				}
			}
			return result;
		}
		SparseMatrix operator-(const SparseMatrix & mat2) const
		{
			SparseMatrix result(*this);
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					result.columnEntries[i][k] -= mat2.columnEntries[i][k];
				}
			}
			return result;
		}
		SparseMatrix & operator=(const SparseMatrix & source)
		{
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					columnEntries[i][k] = source.columnEntries[i][k];
				}
			}
			return *this;
		}
		// The matrix multiply a scalar
		SparseMatrix & operator*=(const T & scalar)
		{
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					columnEntries[i][k] *= scalar;
				}
			}
			return *this;
		}
		SparseMatrix & operator+=(SparseMatrix & mat2)
		{
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					columnEntries[i][k] += mat2.columnEntries[i][k];
				}
			}
			return *this;
		}
		SparseMatrix & operator-=(const SparseMatrix & mat2)
		{
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					columnEntries[i][k] -= mat2.columnEntries[i][k];
				}
			}
			return *this;
		}
		
		// check if the two sparse matrix is the same
		bool operator==(SparseMatrix & mat2)
		{
			if (numRows != mat2.GetNumRows() || numColumns != mat2.GetNumColumns()) {
				return false;
			}
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0, k2; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					if (mat2.CheckEntry(i, j, k2))
						if (columnEntries[i][k] != mat2.columnEntries[i][k2])
							return false;
					else {
						columnEntries[i][k] = 0;
					}
				}
			}
			return true;
		}

		//+ Sparse Matrx Multiply vector or Dense Matrix
		// result = A * vector
		//! 注意矩阵左乘，且维度相同
		void MultiplyVector(const T * vec, T * result) const
		{
			for (size_t i = 0, j; i < numRows; ++i) {
				result[i] = 0;
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					result[i] += columnEntries[i][k] * vec[j];
				}
			}
		}
		// result += A * vector
		void MultiplyVectorAdd(const T * vec, T * result) const
		{
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					result[i] += columnEntries[i][k] * vec[j];
				}
			}
		}
		// result = vector * A
		void TransposeMultiplyVector(const T * vec, T * result) const
		{
			for (size_t j = 0; j < numColumns; ++j) {
				result[j] = 0;
			}
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					result[j] += vec[j] * columnEntries[i][k];
				}
			}
		}
		// result += vector * A
		void TransposeMultiplyVectorAdd(const T * vec, T * result) const
		{
			for (size_t i = 0, j; i < numRows; ++i) {
				for (size_t k = 0; k < rowLength[i]; ++k) {
					j = columnIndices[i][k];
					result[j] += vec[j] * columnEntries[i][k];
				}
			}
		}

		// result = A * denseMatrix
		// A is a (numRows * numColumns) Sparse Matrix 
		// denseMatrix is a (numDenseRows x numDenseColumns) dense matrix
		// result is a (numRows x numDenseColumns) matrix
		//! Notice, numColumns = numDenseRows
		void MultiplyMatrix(size_t numDenseRows, size_t numDenseColumns, const T * denseMatrix, T * result) const
		{
			for (size_t i = 0; i < numRows; ++i) {
				for (size_t j = 0; j < numDenseColumns; j++) {
					result[i*numDenseColumns + j] = 0;
					for (size_t k = 0, t; k < rowLength[i]; ++k) {
						t = columnIndices[i][k];
						//+ result[i][j] += mat[i][t] * dense[t][j]
						//+ mat[i][t] == columnEntries[i][k]
						result[i*numDenseColumns + j] += columnEntries[i][k] * denseMatrix[t*numDenseColumns + j];
					}
				}
			}
		}
		// result += A * denseMatrix
		void MultiplyMatrixAdd(size_t numDenseRows, size_t numDenseColumns, const T * denseMatrix, T * result) const
		{
			for (size_t i = 0; i < numRows; ++i) {
				for (size_t j = 0; j < numDenseColumns; j++) {
					for (size_t k = 0, t; k < rowLength[i]; ++k) {
						t = columnIndices[i][k];
						result[i*numDenseColumns + j] += columnEntries[i][k] * denseMatrix[t*numDenseColumns + j];
					}
				}
			}
		}

		//TODO result = denseMatrix * A & result += denseMatrix * A
		// A 是一个(numRows * numColumns) Sparse Matrix 
		// denseMatrix 是一个 (numDenseRows x numDenseColumns)的矩阵
		// 结果 result 是一个 (numDenseRows x numColumns)的矩阵
		// 注意 numDenseColumns = numRows
		
		// Get Sparse Matrix to array
		void GetEntireSparseMatrixPtr(T * ptr)
		{
			for (size_t i = 0; i < numRows; ++i) {
				size_t j = 0, k = 0, idx;
				while (j < numColumns) {
					idx = i * numColumns + j;
					if (k < rowLength[i] && j == columnIndices[i][k]) {
						ptr[idx] = columnEntries[i][k];
						k += 1;
					}
					else {
						ptr[idx] = 0;
					}
					++j;
				}
			}
		}

		// print matrix
		// print console window in default
		void PrintEntireSparseMatrix(ostream & out = std::cout)
		{
			for (size_t i = 0; i < numRows; ++i) {
				//out << "Row " << i << "#  " << "Length: " << rowLength[i] << endl;
				size_t j = 0, k = 0;
				while (j < numColumns) {
					if (k < rowLength[i] && j == columnIndices[i][k]) {
						out << columnEntries[i][k] << ' ';
						k += 1;
					}
					else {
						out << '0' << ' ';
					}
					++j;
				}
				out << endl;
			}
		}
		
		// overload stream output operator <<
		friend ostream & operator<<(ostream & out, SparseMatrix & source)
		{
			source.PrintEntireSparseMatrix(out);
			return out;
		}

	private:
		// Allocate memory
		void Allocate()
		{
			map_mat.clear();
			map_col_t emptyRow;
			for (size_t i = 0; i < numRows; ++i) {
				map_mat.push_back(emptyRow);
			}
			rowLength.resize(numRows);
			columnIndices.resize(numRows);
			columnEntries.resize(numRows);
		}

	private:
		size_t							numRows;				// 矩阵行数
		size_t							numColumns;				// 矩阵实际列数
		vector< vector<size_t> >		columnIndices;			// 矩阵每一行非零元的位置，由小到大顺序
		vector< vector<T> >				columnEntries;			// 矩阵每一行的非零元
		vector<size_t>					rowLength;				// 矩阵每一行的长度
		map_mat_t						map_mat;				// 矩阵每一行非零元位置和序号映射
	};
}