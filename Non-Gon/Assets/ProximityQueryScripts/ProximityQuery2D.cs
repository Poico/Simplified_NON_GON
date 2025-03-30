using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using static Utility;

public class ProximityQuery2D : MonoBehaviour
{

    //public static bool AABB_AABB2D(Vector3 positionAABB1, Vector3 positionAABB2, double length1, double length2, double height1, double height2)
    //{
    //    Debug.Log("AABB_AABB2D");
    //    bool res = false;
    //    double x1 = positionAABB1.x;
    //    double y1 = positionAABB1.y;
    //    double x2 = positionAABB2.x;
    //    double y2 = positionAABB2.y;
    //    if (x1 < (x2 + length2) && (x1 + length1) > x2 && y1 < (y2 + height2) && (y1 + height1) > y2)
    //    {
    //        res = true;
    //
    //    }
    //    return res;
    //}

    //public static bool OBB_OBB2D(Vector3[] normals1, Vector3[] normals2, Vector3[] corners1, Vector3[] corners2)
    //{
    //    Debug.Log("OBB_OBB2D");
    //    for (int i = 0; i < normals1.Length; i++)
    //    {
    //        float shape1Min = 0, shape1Max = 0, shape2Min = 0, shape2Max = 0;
    //        double[] shape1 = SAT(normals1[i], corners1, shape1Min, shape1Max);
    //        double[] shape2 = SAT(normals1[i], corners2, shape2Min, shape2Max);
    //        if (!Utility.overlaps(shape1[0], shape1[1], shape2[0], shape2[1]))
    //        {
    //            return false;
    //        }
    //    }
    //    for (int i = 0; i < normals2.Length; i++)
    //    {
    //        float shape1Min = 0, shape1Max = 0, shape2Min = 0, shape2Max = 0;
    //        double[] shape1 = SAT(normals2[i], corners1, shape1Min, shape1Max);
    //        double[] shape2 = SAT(normals2[i], corners2, shape2Min, shape2Max);
    //        if (!Utility.overlaps(shape1[0], shape1[1], shape2[0], shape2[1]))
    //        {
    //            return false;
    //        }
    //    }
    //    return true;
    //}

    public static float[] characteristicPolynomial(GameObject Ellipse1, GameObject Ellipse2)
    {
        GeometryCreator variables1 = Ellipse1.GetComponent<GeometryCreator>();
        float xradius1 = variables1.xradius; // Ellipse semi axis
        float yradius1 = variables1.yradius;
        float xEllipse1 = Ellipse1.transform.position.x;
        float yEllipse1 = Ellipse1.transform.position.y;
        GeometryCreator variables2 = Ellipse2.GetComponent<GeometryCreator>();
        float xradius2 = variables2.xradius; // Ellipse semi axis
        float yradius2 = variables2.yradius;
        float xEllipse2 = Ellipse2.transform.position.x;
        float yEllipse2 = Ellipse2.transform.position.y;

        float[,] matrixA = new float[3, 3];
        float[,] matrixB = new float[3, 3];


        //Ellipse 1
        float theta1 = Utility.getRotationZ(Ellipse1) * Mathf.Deg2Rad;
        float sinTheta1 = Mathf.Sin(theta1);
        float cosTheta1 = Mathf.Cos(theta1);
        float aaEllipse1 = xradius1 * xradius1;
        float bbEllipse1 = yradius1 * yradius1;
        float A1 = aaEllipse1 * sinTheta1 * sinTheta1 + bbEllipse1 * cosTheta1 * cosTheta1;
        float B1 = (bbEllipse1 - aaEllipse1) * sinTheta1 * cosTheta1;
        float C1 = aaEllipse1 * cosTheta1 * cosTheta1 + bbEllipse1 * sinTheta1 * sinTheta1;
        float D1 = - A1 * xEllipse1 - B1 * yEllipse1;
        float E1 = - B1 * xEllipse1 - C1 * yEllipse1;
        float F1 = A1 * xEllipse1 * xEllipse1 + 2f * B1 * xEllipse1 * yEllipse1 + C1 * yEllipse1 * yEllipse1 - (aaEllipse1 * bbEllipse1);
        matrixA[0, 0] = A1;
        matrixA[0, 1] = B1;
        matrixA[0, 2] = D1;
        matrixA[1, 0] = B1;
        matrixA[1, 1] = C1;
        matrixA[1, 2] = E1;
        matrixA[2, 0] = D1;
        matrixA[2, 1] = E1;
        matrixA[2, 2] = F1;

        //Debug.Log(A1 + " " + B1 + " " + C1 + " " + D1 + " " + E1 + " " + F1);

        //Ellipse 2
        float theta2 = Utility.getRotationZ(Ellipse2) * Mathf.Deg2Rad;
        float sinTheta2 = Mathf.Sin(theta2);
        float cosTheta2 = Mathf.Cos(theta2);
        float aaEllipse2 = xradius2 * xradius2;
        float bbEllipse2 = yradius2 * yradius2;
        float A2 = aaEllipse2 * sinTheta2 * sinTheta2 + bbEllipse2 * cosTheta2 * cosTheta2;
        float B2 = (bbEllipse2 - aaEllipse2) * sinTheta2 * cosTheta2;
        float C2 = aaEllipse2 * cosTheta2 * cosTheta2 + bbEllipse2 * sinTheta2 * sinTheta2;
        float D2 = - A2 * xEllipse2 - B2 * yEllipse2;
        float E2 = - B2 * xEllipse2 - C2 * yEllipse2;
        float F2 = A2 * xEllipse2 * xEllipse2 + 2f * B2 * xEllipse2 * yEllipse2 + C2 * yEllipse2 * yEllipse2 - (aaEllipse2 * bbEllipse2);
        matrixB[0, 0] = A2;
        matrixB[0, 1] = B2;
        matrixB[0, 2] = D2;
        matrixB[1, 0] = B2;
        matrixB[1, 1] = C2;
        matrixB[1, 2] = E2;
        matrixB[2, 0] = D2;
        matrixB[2, 1] = E2;
        matrixB[2, 2] = F2;

        //Debug.Log(A2 + " " + B2 + " " + C2 + " " + D2 + " " + E2 + " " + F2);

        float a3 = Utility.calDetMatrix3x3(matrixA);
        float a2 = (1 / 3f) * (A1 * C1 * F2 + A1 * C2 * F1 + A2 * C1 * F1 + B1 * E1 * D2 + B1 * E2 * D1 + B2 * E1 * D1 + D1 * B1 * E2 + D1 * B2 * E1 + D2 * B1 * E1 -
            D1 * C1 * D2 - D1 * C2 * D1 - D2 * C1 * D1 - B1 * B1 * F2 - B1 * B2 * F1 - B2 * B1 * F1 - A1 * E1 * E2 - A1 * E2 * E1 - A2 * E1 * E1);
        float a1 = (1 / 3f) * (A1 * C2 * F2 + A2 * C1 * F2 + A2 * C2 * F1 + B1 * E2 * D2 + B2 * E1 * D2 + B2 * E2 * D1 + D1 * B2 * E2 + D2 * B1 * E2 + D2 * B2 * E1 -
            D1 * C2 * D2 - D2 * C1 * D2 - D2 * C2 * D1 - B1 * B2 * F2 - B2 * B1 * F2 - B2 * B2 * F1 - A1 * E2 * E2 - A2 * E1 * E2 - A2 * E2 * E1);
        float a0 = Utility.calDetMatrix3x3(matrixB);

        /*
        // This also calculates the characteristic polynomial 
        float[,] m1 = new float[3, 3] { { A2, B1, D1 }, { B2, C1, E1 }, { D2, E1, F1 } };
        float[,] m2 = new float[3, 3] { { A1, B2, D1 }, { B1, C2, E1 }, { D1, E2, F1 } };
        float[,] m3 = new float[3, 3] { { A1, B1, D2 }, { B1, C1, E2 }, { D1, E1, F2 } };
        float[,] m4 = new float[3, 3] { { A1, B2, D2 }, { B1, C2, E2 }, { D1, E2, F2 } };
        float[,] m5 = new float[3, 3] { { A2, B1, D2 }, { B2, C1, E2 }, { D2, E1, F2 } };
        float[,] m6 = new float[3, 3] { { A2, B2, D1 }, { B2, C2, E1 }, { D2, E2, F1 } };

        a2 = (1 / 3f) * (Utility.calDetMatrix3x3(m1) + Utility.calDetMatrix3x3(m2) + Utility.calDetMatrix3x3(m3));
        a1 = (1 / 3f) * (Utility.calDetMatrix3x3(m4) + Utility.calDetMatrix3x3(m5) + Utility.calDetMatrix3x3(m6));
        */

        /*
        Debug.Log("a3 = " + a3);
        Debug.Log("a2 = " + a2);
        Debug.Log("a1 = " + a1);
        Debug.Log("a0 = " + a0);
        Debug.Log("--------------------------------");
        */

        return new float[4] { a3, a2, a1, a0 };
    }

    public static bool Ellipse_Ellipse_Caravantes(GameObject Ellipse1, GameObject Ellipse2)
    {
        float[] characteristicPolynomialValues = characteristicPolynomial(Ellipse1, Ellipse2);

        float a3 = characteristicPolynomialValues[0];
        float a2 = characteristicPolynomialValues[1];
        float a1 = characteristicPolynomialValues[2];
        float a0 = characteristicPolynomialValues[3];

        float s0 = a3 * ((27f * a0 * a0 * a3 * a3) - (18f * a0 * a1 * a2 * a3) + (4f * a0 * a2 * a2 * a2) + (4f * a1 * a1 * a1 * a3) - (a1 * a1 * a2 * a2));

        if (DescartesLawOfSigns.DescartesLawOfSignsThirdDegreePolynomial(a3, a2, a1, a0)) {
            if (s0 >= 0f)
            {
                return false;
            } else
            {
                return true;
            }
        } else
        {
            return true;
        }
    }

    //public static bool Ellipse_Ellipse_Alberich(GameObject Ellipse1, GameObject Ellipse2)
    //{
    //    Debug.Log("Ellipse_Ellipse_Alberich");
    //    float[] characteristicPolynomialValues = characteristicPolynomial(Ellipse1, Ellipse2);
    //
    //    float l3 = characteristicPolynomialValues[0];
    //    float l2 = characteristicPolynomialValues[1];
    //    float l1 = characteristicPolynomialValues[2];
    //    float l0 = characteristicPolynomialValues[3];
    //
    //    float delta1 = l3 * l1 - l2 * l2;
    //    float delta2 = l3 * l0 - l2 * l1;
    //    float delta3 = l2 * l0 - l1 * l1;
    //
    //    float discriminant_P = 4f * delta1 * delta3 - delta2 * delta2;
    //
    //    if (discriminant_P >= 0f)
    //    {
    //        if (l1 > 0f || l2 > 0f)
    //        {
    //            if (discriminant_P > 0f)
    //            {
    //                return false;
    //            } else
    //            {
    //                return true;
    //            }
    //        } else
    //        {
    //            return true;
    //        }
    //    } else
    //    {
    //        return true;
    //    }
    //}
}