/*
This is an implementation of the GSW scheme using the NTL and GMP C++ libraries.

Sources:
NTL C++ library: https://libntl.org/
GMP C++ library: https://gmplib.org/

Articles and references:
GSW scheme: https://eprint.iacr.org/2013/340.pdf

other related articles:
https://web.eecs.umich.edu/~cpeikert/pubs/polyboot.pdf
https://eprint.iacr.org/2021/691.pdf
https://eprint.iacr.org/2020/086.pdf
*/

#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>

using namespace std;
using namespace NTL;

namespace GSW {
    //////////////////////////////////////////////////////////////////////////////
    // forward declaration:
    Vec<ZZ_p> G_inverse(ZZ_p num, long bit_length);

    //////////////////////////////////////////////////////////////////////////////
    void print_VecOfVecs(Vec<Vec<ZZ>> vecOfVecs) {
        for (size_t i = 0; i < vecOfVecs.length(); i++) {
            cout << vecOfVecs[i] << endl;
        }
    }

    void print_Matrix(Mat<ZZ_p> matrix) {
        for (size_t i = 0; i < matrix.NumRows(); i++) {
            for (size_t j = 0; j < matrix.NumCols(); j++) {
                cout << matrix[i][j] << endl;
            }
        }
    }

    Mat<ZZ_p> G_inverse_of_ciphertext(Mat<ZZ_p> ciphertext, long n, long bit_length) {
        // this function  gets a ciphertext as input, which is a matrix of dimension (2 ⨯ 2*bit_length).
        // The function returns the bitdecomposition of the ciphertext.
        // the output is a matrix of dimension (2*bit_length ⨯ 2*bit_length) with binary entries.

        Mat<ZZ_p> bitdecomposed_ciphertext;
        bitdecomposed_ciphertext.SetDims(n * bit_length, n * bit_length);

        Vec<ZZ_p> num_bitdecomposed;
        num_bitdecomposed.SetLength(bit_length);

        ZZ_p temp_num;
        for (size_t i = 0; i < n * bit_length; i++) {
            for (size_t j = 0; j < n; j++) {
                temp_num = ciphertext[j][i];
                num_bitdecomposed = G_inverse(temp_num, bit_length);
                for (size_t k = 0; k < bit_length; k++) {
                    bitdecomposed_ciphertext[j * bit_length + k][i] = num_bitdecomposed[k];
                }
            }
        }

        return bitdecomposed_ciphertext;
    }

    Vec<ZZ_p> G_inverse(ZZ_p num, long bit_length) {
        // this function gets an integer number and bit-decomposes the number the way that G^-1 works:

        Vec<ZZ_p> num_bitdecomposed;
        num_bitdecomposed.SetLength(bit_length);

        ZZ num_in_ZZ;
        for (size_t i = 0; i < bit_length; i++) {
            conv(num_in_ZZ, num);
            num_bitdecomposed[i] = bit(num_in_ZZ, i);
        }

        return num_bitdecomposed;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Implementation of the encryption algorithm for the GSW scheme
    Mat<ZZ_p> encrypt(int message, Mat<ZZ_p> A, Mat<ZZ_p> R, long n, long bit_length) {

        long l = bit_length;

        ZZ_p m;
        m = ZZ_p(message);

        // calculating A * R:
        Mat<ZZ_p> A_milttipliedBy_R;
        A_milttipliedBy_R.SetDims(n, n * l);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n * l; j++) {
                A_milttipliedBy_R[i][j] = ZZ_p(0);
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n * l; j++) {
                for (int k = 0; k < n * l; k++) {
                    A_milttipliedBy_R[i][j] += A[i][k] * R[k][j];
                }
            }
        }

        // creating matrix G:
        Mat<ZZ_p> G;
        G.SetDims(n, n * l);
        ZZ_p temp = ZZ_p(1);

        //creating the matrix G:
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n * l; j++) {
                G[i][j] = ZZ_p(0);
            }
        }

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < l; j++) {
                G[i][i * l + j] = temp;
                mul(temp, temp, ZZ_p(2));
            }
            temp = ZZ_p(1);
        }

        // calculating m . G:
        Mat<ZZ_p> m_milttipliedBy_G;
        m_milttipliedBy_G.SetDims(n, n * l);

        if (message == 1)
            m_milttipliedBy_G = G;
        else { // in this case (i.e., if m == 0): m . G = 0.
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2 * l; j++) {
                    m_milttipliedBy_G[i][j] = ZZ_p(0);
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////
        // Encrypting a message (where m can be 0 or 1) using the GSW scheme:
        // i.e., ciphertext = A . R + m . G:
        Mat<ZZ_p> ciphertext;
        ciphertext.SetDims(n, n * l);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n * l; j++) {
                ciphertext[i][j] = A_milttipliedBy_R[i][j] + m_milttipliedBy_G[i][j];
            }
        }

        return ciphertext;
    }

    ZZ_p dot_product(Vec<ZZ_p> v1, Vec<ZZ_p> v2) {

        if (v1.length() != v2.length()) {
            cout << "The lengths of vectors are different!" << endl;
            cout << "v1.length(): " << v1.length() << "          v2.length(): " << v2.length() << endl;
        }

        ZZ_p result = ZZ_p(0);
        ZZ_p temp;
        for (size_t i = 0; i < v1.length(); i++) {
            temp = v1[i] * v2[i];
            result += temp;
        }

        return result;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Implementation of the decryption algorithm for the RGSW scheme
    int decrypt(Vec<ZZ_p> sk, Mat<ZZ_p> ciphertext, long n, long bit_length, ZZ modulous) {
        int decrypted_message;

        decrypted_message = 0; // initializing the variabele.
        Vec<ZZ_p> penultimate_column_of_c, last_column_of_c;
        Vec<ZZ_p> first_column_of_c, second_column_of_c;

        for (size_t i = 0; i < n; i++) {
            first_column_of_c.append(ciphertext[i][0]);
            second_column_of_c.append(ciphertext[i][1]);

            penultimate_column_of_c.append(ciphertext[i][2 * bit_length - 2]);
            last_column_of_c.append(ciphertext[i][2 * bit_length - 1]);
        }

        ZZ_p result1, result2, result3, result4;
        ZZ result3zz;
        long comparison_result;

        result1 = dot_product(sk, first_column_of_c);
        result2 = dot_product(sk, second_column_of_c);
        result3 = dot_product(sk, penultimate_column_of_c);
        result4 = dot_product(sk, last_column_of_c);

        ZZ p_dividedBy_2;;
        div(p_dividedBy_2, modulous, ZZ(2));
        conv(result3zz, result3);
        comparison_result = compare(result3zz, p_dividedBy_2);

        if (comparison_result > 0)
            decrypted_message = 1;
        else
            decrypted_message = 0;

        return decrypted_message;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Implementation of the addition gate for the GSW scheme
    Mat<ZZ_p> add(Mat<ZZ_p> c1, Mat<ZZ_p> c2) {
        Mat<ZZ_p> addtion_result;

        if ((c1.NumRows() != c2.NumRows()) || (c1.NumCols() != c2.NumCols())) {
            cout << "Error in the add function: the dimensions of the ciphertext do not matach!" << endl;
        }

        addtion_result.SetDims(c1.NumRows(), c1.NumCols());

        for (size_t i = 0; i < c1.NumRows(); i++) {
            for (size_t j = 0; j < c1.NumCols(); j++) {
                addtion_result[i][j] = c1[i][j] + c2[i][j];
            }
        }

        return addtion_result;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Implementation of the multiplication gate for the GSW scheme
    Mat<ZZ_p> multiply(Mat<ZZ_p> c1, Mat<ZZ_p> c2, long n, long bit_length) {

        if ((c1.NumRows() != c2.NumRows()) || (c1.NumCols() != c2.NumCols())) {
            cout << "Error in the multiply function: the dimensions of the ciphertext do not matach!" << endl;
        }

        Mat<ZZ_p> multiplication_result;
        multiplication_result.SetDims(c1.NumRows(), c2.NumCols());

        // calculating G_inverse (c2):
        Mat<ZZ_p> G_inverseOf_c2 = G_inverse_of_ciphertext(c2, n, bit_length);

        // initializing all elements of multiplication_result to zero.
        multiplication_result.SetDims(c1.NumRows(), c2.NumCols());
        for (size_t i = 0; i < c1.NumRows(); i++) {
            for (size_t j = 0; j < c2.NumCols(); j++) {
                multiplication_result[i][j] = ZZ_p(0);
            }
        }

        // calculating c1 . G_inverse(c2):
        ZZ_p temp;
        for (size_t i = 0; i < c1.NumRows(); i++) {
            for (size_t j = 0; j < c1.NumCols(); j++) {
                for (size_t k = 0; k < G_inverseOf_c2.NumRows(); k++) {
                    //cout << "i: " << i << ", j: " << j << ", k: " << k << endl;
                    temp = c1[i][k] * G_inverseOf_c2[k][j];
                    multiplication_result[i][j] += temp;
                }
            }
        }

        return multiplication_result;
    }


    void myGSW() {

        ZZ p;
        long bit_length = 128;
        GenPrime(p, bit_length);
        cout << "The prime modulous (i.e., p) is: " << p;
        cout << "   (it is a " << bit_length << " bits prime number)" << endl;
        ZZ_p::init(p);

        long n = 4; // matrix A in the GSW scheme has dimensions n ⨯ n*l. This 'n' variable is used for showing the number of rows in the matrix A.
        cout << "Please wait......., some computations might take some time...";

        //////////////////////////////////////////////////////////////////////////////
        // generating the secret key:
        Vec<ZZ_p> s, sk;
        s.SetLength(n - 1);
        sk.SetLength(n);

        for (size_t i = 0; i < s.length() - 1; i++) {
            random(s[i]);
            sk[i] = s[i];
        }
        sk[n - 1] = ZZ_p(1);

        // generating an error vector:
        Vec<ZZ_p> error_vec;
        error_vec.SetLength(n * bit_length);
        // setting the modulous to 2, for generating a small error vector:
        ZZ_p::init(ZZ(2));
        for (size_t i = 0; i < n * bit_length; i++) {
            random(error_vec[i]);
        }

        // setting the modulous back to p:
        ZZ_p::init(p);

        Mat<ZZ_p> A_temp, A;
        A_temp.SetDims(n - 1, n * bit_length);
        A.SetDims(n, n * bit_length);

        ZZ_p random_num;

        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n * bit_length; j++) {
                random(random_num);
                A_temp[i][j] = random_num;
                A[i][j] = random_num;
            }
        }

        Vec<ZZ_p> s_multipliedBy_A_temp;
        s_multipliedBy_A_temp.SetLength(A_temp.NumCols());
        for (size_t i = 0; i < A_temp.NumCols(); i++) {
            s_multipliedBy_A_temp[i] = ZZ_p(0);
        }

        ZZ_p mult_temp;
        for (size_t i = 0; i < s.length(); i++) {
            for (size_t j = 0; j < A_temp.NumCols(); j++) {
                for (size_t k = 0; k < A_temp.NumRows(); k++) {
                    mult_temp = s[k] * A_temp[k][j];
                    s_multipliedBy_A_temp[j] += mult_temp;
                }
            }
        }

        // creating the public key, i.e., matrix A,
        // which is calculated as: e - s*A_temp:
        for (int i = n - 1; i < n; i++) {
            for (int j = 0; j < n * bit_length; j++) {
                A[i][j] = error_vec[j] - s_multipliedBy_A_temp[j];
            }
        }

        // constructing the random matrix R (as defined in the GSW scheme) with binary entries:
        // setting the modulous to 2, for generating a small random numbers:
        ZZ_p::init(ZZ(2));
        ZZ_p::init(ZZ(2));
        Mat<ZZ_p> R;
        R.SetDims(n * bit_length, n * bit_length);
        ZZ_p temp_num;
        for (size_t i = 0; i < n * bit_length; i++) {
            for (size_t j = 0; j < n * bit_length; j++) {
                random(temp_num);
                R[i][j] = temp_num;
            }
        }

        // setting the modulous back to p:
        ZZ_p::init(p);


        //////////////////////////////////////////////////////////////////////////////
        // creating two ciphertexts and testing the add and multiply operations on ciphertexts:
        int message1 = 1, message2 = 0;
        Mat<ZZ_p> ciphertext1, ciphertext2;
        Mat<ZZ_p> addition_result_ct, multiplication_result_ct;


        ciphertext1 = encrypt(message1, A, R, n, bit_length);

        ciphertext2 = encrypt(message2, A, R, n, bit_length);

        //cout << "ciphertext1: " << endl;
        //cout << ciphertext1 << endl;

        //cout << "ciphertext2: " << endl;
        //cout << ciphertext2 << endl;

        addition_result_ct = add(ciphertext1, ciphertext2);

        multiplication_result_ct = multiply(ciphertext1, ciphertext2, n, bit_length);

        //cout << "\naddition_result: " << endl;
        //cout << addition_result << endl;

        //cout << "multiplication_result: " << endl;
        //cout << multiplication_result << endl;

        //////////////////////////////////////////////////////////////////////////////
        // decrypting the result:
        int decryption_result1, decryption_result2;

        decryption_result1 = decrypt(sk, ciphertext1, n, bit_length, p);
        decryption_result2 = decrypt(sk, ciphertext2, n, bit_length, p);

        cout << "\n\nmessage 1 is: " << message1 << endl;
        cout << "message 2 is: " << message2 << endl;
        cout << "decryption_result1 is: " << decryption_result1 << endl;
        cout << "decryption_result2 is: " << decryption_result2 << endl;

    }
}


int main() {
    
    cout << "################################################################" << endl;
    cout << "#------------- Welcome to GSW Implementation! -------------#" << endl;
    cout << "################################################################" << endl;

    GSW::myGSW();

    cout << endl;
    cout << "################################################################" << endl;
    cout << "#----------------------- End of Program -----------------------#" << endl;
    cout << "################################################################" << endl;

    return 0;
}