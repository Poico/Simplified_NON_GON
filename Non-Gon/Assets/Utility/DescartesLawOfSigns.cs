using System.Collections.Generic;
using UnityEngine;

public class DescartesLawOfSigns : MonoBehaviour
{
    // Check for 2 change of signs according to Descartes Law of Signs
    public static bool DescartesLawOfSignsFourthDegreePolynomial(float a4, float a3, float a2, float a1, float a0)
    {
        // Patterns based on signs of coefficients
        // Examples for patterns (for illustration purposes)
        //  -, +, +, +, -
        //  -, -, +, -, -
        //  -, -, +, +, -
        //  -, -, -, +, -
        //  -, +, -, -, -
        //  -, +, +, -, -
        //  -, +, 0, -, -
        //  -, 0, +, -, -
        //  -, +, -, 0, -
        //  -, -, 0, +, -
        //  -, -, +, 0, -
        //  -, 0, -, +, -
        //  -, +, 0, 0, -
        //  -, 0, +, 0, -
        //  -, 0, 0, +, -
        //  -, +, 0, +, -
        //  -, +, +, 0, -
        //  -, 0, +, +, -

        bool pattern1to18 = (a4 < 0 && a3 > 0 && a2 > 0 && a1 > 0 && a0 < 0) || (a4 < 0 && a3 < 0 && a2 > 0 && a1 < 0 && a0 < 0) || (a4 < 0 && a3 < 0 && a2 > 0 && a1 > 0 && a0 < 0) ||
            (a4 < 0 && a3 < 0 && a2 < 0 && a1 > 0 && a0 < 0) || (a4 < 0 && a3 > 0 && a2 < 0 && a1 < 0 && a0 < 0) || (a4 < 0 && a3 > 0 && a2 > 0 && a1 < 0 && a0 < 0) ||
            (a4 < 0 && a3 > 0 && a2 == 0 && a1 < 0 && a0 < 0) || (a4 < 0 && a3 == 0 && a2 > 0 && a1 < 0 && a0 < 0) || (a4 < 0 && a3 > 0 && a2 < 0 && a1 == 0 && a0 < 0) ||
            (a4 < 0 && a3 < 0 && a2 == 0 && a1 > 0 && a0 < 0) || (a4 < 0 && a3 < 0 && a2 > 0 && a1 == 0 && a0 < 0) || (a4 < 0 && a3 == 0 && a2 < 0 && a1 > 0 && a0 < 0) ||
            (a4 < 0 && a3 > 0 && a2 == 0 && a1 == 0 && a0 < 0) || (a4 < 0 && a3 == 0 && a2 > 0 && a1 == 0 && a0 < 0) || (a4 < 0 && a3 == 0 && a2 == 0 && a1 > 0 && a0 < 0) ||
            (a4 < 0 && a3 > 0 && a2 == 0 && a1 > 0 && a0 < 0) || (a4 < 0 && a3 > 0 && a2 > 0 && a1 == 0 && a0 < 0) || (a4 < 0 && a3 == 0 && a2 > 0 && a1 > 0 && a0 < 0);

        return pattern1to18;
    }

    // Check for 2 change of signs according to Descartes Law of Signs
    public static bool DescartesLawOfSignsThirdDegreePolynomial(float a3, float a2, float a1, float a0)
    {
        // Patterns based on signs of coefficients
        // Examples for patterns (for illustration purposes)
        // -, 0, +, -
        // -, +, 0, -
        // -, +, -, -
        // -, +, +, -
        // -, -, +, -

        return (a3 < 0 && a2 == 0 && a1 > 0 && a0 < 0) || (a3 < 0 && a2 > 0 && a1 == 0 && a0 < 0) || (a3 < 0 && a2 > 0 && a1 < 0 && a0 < 0) ||
            (a3 < 0 && a2 > 0 && a1 > 0 && a0 < 0) || (a3 < 0 && a2 < 0 && a1 > 0 && a0 < 0);
    }
}