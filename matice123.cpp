// matice.cpp: Definuje vstupní bod pro aplikaci.
//


#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
using namespace std;

class Matice {
    int radky, sloupce;
    double mat[5][5];

public:
    Matice(int r, int s) : radky(r), sloupce(s) {}

    void nastavHodnoty() {
        cout << "Zadejte hodnoty pro matici (" << radky << "x" << sloupce << "):\n";
        for (int i = 0; i < radky; i++) {
            for (int j = 0; j < sloupce; j++) {
                cout << "Prvek [" << i + 1 << "][" << j + 1 << "]: ";
                cin >> mat[i][j];

                if (cin.fail()) {
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                    cout << "Neplatny vstup! Zadejte cislo.\n";
                    j--;
                }
            }
        }
    }

    void vypisMatici() const {
        for (int i = 0; i < radky; i++) {
            for (int j = 0; j < sloupce; j++) {
                cout << mat[i][j] << " ";
            }
            cout << endl;
        }
    }

    Matice soucet(const Matice& B) const {
        Matice vysledek(radky, sloupce);
        for (int i = 0; i < radky; i++) {
            for (int j = 0; j < sloupce; j++) {
                vysledek.mat[i][j] = mat[i][j] + B.mat[i][j];
            }
        }
        return vysledek;
    }

    Matice rozdil(const Matice& B) const {
        Matice vysledek(radky, sloupce);
        for (int i = 0; i < radky; i++) {
            for (int j = 0; j < sloupce; j++) {
                vysledek.mat[i][j] = mat[i][j] - B.mat[i][j];
            }
        }
        return vysledek;
    }

    Matice nasobSkalarem(double skalar) const {
        Matice vysledek(radky, sloupce);
        for (int i = 0; i < radky; i++) {
            for (int j = 0; j < sloupce; j++) {
                vysledek.mat[i][j] = mat[i][j] * skalar;
            }
        }
        return vysledek;
    }

    Matice nasobitMatici(const Matice& B) const {
        if (sloupce != B.radky) {
            throw runtime_error("Matice nelze nasobit, nesouhlasi rozmery.");
        }
        Matice vysledek(radky, B.sloupce);
        for (int i = 0; i < radky; i++) {
            for (int j = 0; j < B.sloupce; j++) {
                vysledek.mat[i][j] = 0;
                for (int k = 0; k < sloupce; k++) {
                    vysledek.mat[i][j] += mat[i][k] * B.mat[k][j];
                }
            }
        }
        return vysledek;
    }

    Matice transponovana() const {
        Matice vysledek(sloupce, radky);
        for (int i = 0; i < radky; i++) {
            for (int j = 0; j < sloupce; j++) {
                vysledek.mat[j][i] = mat[i][j];
            }
        }
        return vysledek;
    }

    double determinant() const {
        if (radky != sloupce) {
            throw runtime_error("Determinant lze pocitat jen pro ctvercove matice.");
        }
        if (radky == 2) {
            return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
        }
        if (radky == 3) {
            return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
        }
        if (radky == 4) {
            double det = 0;
            for (int i = 0; i < 4; i++) {
                Matice subMat(3, 3);
                for (int j = 1; j < 4; j++) {
                    int colIndex = 0;
                    for (int k = 0; k < 4; k++) {
                        if (k == i) continue;
                        subMat.mat[j - 1][colIndex] = mat[j][k];
                        colIndex++;
                    }
                }
                double subDet = subMat.determinant();
                det += (i % 2 == 0 ? 1 : -1) * mat[0][i] * subDet;
            }
            return det;
        }
        throw runtime_error("Determinant je podporovan jen do velikosti 4x4.");
    }
    
        int hodnost() const {
        Matice A = *this;
        int r = radky;
        int c = sloupce;
        int rank = std::min(r, c);

        for (int i = 0; i < rank; i++) {
            if (A.mat[i][i] == 0) {
                bool swapFound = false;
                for (int j = i + 1; j < r; j++) {
                    if (A.mat[j][i] != 0) {
                        swapFound = true;

                        for (int k = 0; k < c; k++) {
                            double temp = A.mat[i][k];
                            A.mat[i][k] = A.mat[j][k];
                            A.mat[j][k] = temp;
                        }
                        break;
                    }
                }
                if (!swapFound) {
                    rank--;
                }
            }

            if (A.mat[i][i] != 0) {
                for (int j = i + 1; j < r; j++) {
                    double ratio = A.mat[j][i] / A.mat[i][i];
                    if (A.mat[i][i] != 0) {
                        for (int k = i; k < c; k++) {
                            A.mat[j][k] -= ratio * A.mat[i][k];
                        }
                    }
                }
            }
        }
        return rank;
    }

    Matice inverzniMatice() const {
        if (radky != sloupce) {
            throw runtime_error("Inverzni matice lze vypocitat pouze pro ctvercove matice.");
        }
        Matice A = *this;
        Matice B(radky, sloupce);


        for (int i = 0; i < radky; i++) {
            for (int j = 0; j < sloupce; j++) {
                B.mat[i][j] = (i == j) ? 1 : 0;
            }
        }

        for (int i = 0; i < radky; i++) {
            double diagElem = A.mat[i][i];
            if (diagElem == 0) {
                throw runtime_error("Matice nema inverzi.");
            }
            for (int j = 0; j < sloupce; j++) {
                A.mat[i][j] /= diagElem;
                B.mat[i][j] /= diagElem;
            }
            for (int k = 0; k < radky; k++) {
                if (k != i) {
                    double factor = A.mat[k][i];
                    for (int j = 0; j < sloupce; j++) {
                        A.mat[k][j] -= factor * A.mat[i][j];
                        B.mat[k][j] -= factor * B.mat[i][j];
                    }
                }
            }
        }
        return B;
    }
};
int main() {
    int volba;
    int r, s;
    cout << "Zadejte pocet radku matice (max 5): ";
    cin >> r;
    cout << "Zadejte pocet sloupcu matice (max 5): ";
    cin >> s;

    Matice A(r, s);
    A.nastavHodnoty();

    while (true) {
        cout << "\nVyberte operaci:\n";
        cout << "1. Soucet matic\n";
        cout << "2. Rozdil matic\n";
        cout << "3. Nasobeni skalarem\n";
        cout << "4. Nasobeni matic\n";
        cout << "5. Transponovana matice\n";
        cout << "6. Determinant (pro matice 2x2 az 4x4)\n";
        cout << "7. Hodnost matice\n";
        cout << "8. Inverzni matice\n";
        cout << "9. Konec\n";
        cin >> volba;

        if (volba == 9) {
            break;
        }

        switch (volba) {
        case 1: {
            Matice B(r, s);
            B.nastavHodnoty();
            Matice C = A.soucet(B);
            C.vypisMatici();
            break;
        }
        case 2: {
            Matice B(r, s);
            B.nastavHodnoty();
            Matice C = A.rozdil(B);
            C.vypisMatici();
            break;
        }
        case 3: {
            double skalar;
            cout << "Zadejte skalar: ";
            cin >> skalar;
            Matice C = A.nasobSkalarem(skalar);
            C.vypisMatici();
            break;
        }
        case 4: {
            Matice B(s, r);
            B.nastavHodnoty();
            try {
                Matice C = A.nasobitMatici(B);
                C.vypisMatici();
            }
            catch (runtime_error& e) {
                cout << e.what() << endl;
            }
            break;
        }
        case 5: {
            Matice C = A.transponovana();
            C.vypisMatici();
            break;
        }
        case 6: {
            try {
                double det = A.determinant();
                cout << "Determinant: " << det << endl;
            }
            catch (runtime_error& e) {
                cout << e.what() << endl;
            }
            break;
        }
        case 7: {
            int rank = A.hodnost();
            cout << "Hodnost matice: " << rank << endl;
            break;
        }
        case 8: {
            try {
                Matice C = A.inverzniMatice();
                cout << "Inverzni matice:\n";
                C.vypisMatici();
            }
            catch (runtime_error& e) {
                cout << e.what() << endl;
            }
            break;
        }
        default:
            cout << "Neplatna volba.\n";
        }
    }

    return 0;
}