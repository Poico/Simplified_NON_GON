using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using Complex = System.Numerics.Complex;

public class Utility : MonoBehaviour
{ 

    public static float CubeRoot(float d)
    {
    if (d < 0.0f) {
        return -Mathf.Pow(-d, 1f / 3f);
    }
    else {
        return Mathf.Pow(d, 1f / 3f);
    }
    }

    public static double CubeRootd(double d)
    {
    if (d < 0.0d) {
        return -Math.Pow(-d, 1d / 3d);
    }
    else {
        return Math.Pow(d, 1d / 3d);
    }
    }

    public static Complex CubeRootc (Complex c)
    {
        float phi = (float)c.Phase;
        float r = (float)c.Magnitude;
        phi /= 3.0f;
        r = CubeRoot(r);
        return r*Mathf.Cos(phi) + Complex.ImaginaryOne*r*Mathf.Sin(phi);
    }

    public static float[] quartic_roots (float a, float b, float c, float d, float e)
    {
        double s0 = c*c - 3.0*b*d + 12.0*a*e;
        double s1 = 2.0*c*c*c - 9.0*b*c*d + 27.0*b*b*e + 27.0*a*d*d - 72.0*a*c*e;
        double p = (8.0*a*c - 3.0*b*b)/(8.0*a*a);
        double q = (b*b*b - 4.0*a*b*c + 8.0*a*a*d)/(8.0*a*a*a);

        double Q = CubeRootd(0.5*(s1 + Math.Sqrt(s1*s1 - 4.0*s0*s0*s0)));
        double S;
        if(double.IsNaN(Q))
        {
            double phi = Math.Acos(s1/(2.0*Math.Sqrt(s0*s0*s0)));
            S = 0.5*Math.Sqrt(-2.0/3.0*p + 2.0/(3.0*a)*Math.Sqrt(s0)*Math.Cos(phi/3.0));
        }
        else
        {
            S = 0.5*Math.Sqrt(-(2.0/3.0)*p + 1.0/(3.0*a)*(Q + s0/Q));
        }

        float[] sols = new float[4];

        sols[0] = System.Convert.ToSingle(-0.25*b/a - S + 0.5f*Math.Sqrt(-4.0*S*S - 2.0*p + q/S));
        sols[1] = System.Convert.ToSingle(-0.25*b/a - S - 0.5f*Math.Sqrt(-4.0*S*S - 2.0*p + q/S));
        sols[2] = System.Convert.ToSingle(-0.25*b/a + S + 0.5f*Math.Sqrt(-4.0*S*S - 2.0*p - q/S));
        sols[3] = System.Convert.ToSingle(-0.25*b/a + S - 0.5f*Math.Sqrt(-4.0*S*S - 2.0*p - q/S));
        
        return sols;
    }

    public static float[] quartic_roots_ (float a, float b, float c, float d, float e)
    {
        float[] roots = new float[4];

        double z = b*b - 3.0*a*c +12.0*d;
        double z0 = 2.0*b*b*b - 9.0*a*b*c + 27.0*c*c + 27.0*a*a*d - 72.0*b*d;
        double z2 = z0 + Math.Sqrt(-4.0*z*z*z + z0*z0);
        double z1 = Math.Pow(2.0, 1.0/3.0)*z/(3.0*Math.Pow(z2,1.0/3.0));
        double z3 = Math.Sqrt(a*a/4.0 - 2.0*b/3.0 + z1 + Math.Pow(z2/54.0,1.0/3.0));
        double z4 = -a*a*a + 4.0*a*b - 8.0*c;
        double z5 = 0.5*a*a - 4.0*b/3.0 - z1 - Math.Pow(z2/54.0,1.0/3.0);
        
        roots[0] = System.Convert.ToSingle(-a/4.0 - 0.5*z3 - 0.5*Math.Sqrt(z5 - z4/(4.0*z3)));
        roots[1] = System.Convert.ToSingle(-a/4.0 - 0.5*z3 + 0.5*Math.Sqrt(z5 - z4/(4.0*z3)));
        roots[2] = System.Convert.ToSingle(-a/4.0 + 0.5*z3 - 0.5*Math.Sqrt(z5 + z4/(4.0*z3)));
        roots[3] = System.Convert.ToSingle(-a/4.0 + 0.5*z3 + 0.5*Math.Sqrt(z5 + z4/(4.0*z3)));

        return roots;
    }

    public static float Xroot(float a, float x)
    {
        float i = 1f;
        if (a < 0f)
            i = -1f;
        return (i * Mathf.Exp(Mathf.Log(a * i) / x));
    }

    public static Complex[] SolveQuarticEquation(float a4, float a3, float a2, float a1, float a0)
    {
        float[] coeffs = { a0 / a4, a1 / a4, a2 / a4, a3 / a4, 1 }; // Normalize coefficients

        float A = coeffs[3];
        float B = coeffs[2];
        float C = coeffs[1];
        float D = coeffs[0];

        // Calculate coefficients for resolvent cubic equation
        float p = C - (3 * B * B) / 8;
        float q = D + (B * B * B) / 8 - (B * C) / 2;
        float r = -(3 * B * B * B * B) / 256 + (C * B * B) / 16 - (B * D) / 4 + A;

        // Solve resolvent cubic equation
        Complex[] cubicRoots = SolveCubicEquation(1, p, -q, -r);

        float z1 = (float)cubicRoots[0].Real;
        float z2 = (float)cubicRoots[1].Real;
        float z3 = (float)cubicRoots[2].Real;

        Complex[] roots = new Complex[4];

        float sqrt2 = (float)Mathf.Sqrt(2);

        roots[0] = new Complex(-B / 4 + sqrt2 * (float)Math.Sqrt(z1) + 0.5f * (float)Math.Sqrt(-(2 * z1 + p - (B * B) / 4) + 2 * sqrt2 * (z2 + z3)), 0);
        roots[1] = new Complex(-B / 4 - sqrt2 * (float)Math.Sqrt(z1) + 0.5f * (float)Math.Sqrt(-(2 * z1 + p - (B * B) / 4) - 2 * sqrt2 * (z2 + z3)), 0);
        roots[2] = new Complex(-B / 4 + sqrt2 * (float)Math.Sqrt(z2) + 0.5f * (float)Math.Sqrt(-(2 * z2 + p - (B * B) / 4) + 2 * sqrt2 * (z1 + z3)), 0);
        roots[3] = new Complex(-B / 4 - sqrt2 * (float)Math.Sqrt(z2) + 0.5f * (float)Math.Sqrt(-(2 * z2 + p - (B * B) / 4) - 2 * sqrt2 * (z1 + z3)), 0);

        return roots;
    }

    private static Complex[] SolveCubicEquation(float a, float b, float c, float d)
    {
        Complex[] roots = new Complex[3];

        Complex p = (3 * a * c - b * b) / (3 * a * a);
        Complex q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);

        Complex discriminant = q * q / 4 + p * p * p / 27;
        Complex omega = Complex.Exp(new Complex(0, 2 * Math.PI / 3));

        if (discriminant == 0 && q == 0)
        {
            roots[0] = new Complex(-b / (3 * a), 0);
            roots[1] = new Complex(-b / (3 * a), 0);
            roots[2] = new Complex(-b / (3 * a), 0);
        }
        else if (discriminant.Real > 0)
        {
            Complex u = Complex.Pow(-q / 2 + Complex.Sqrt(discriminant), 1.0f / 3.0f);
            Complex v = Complex.Pow(-q / 2 - Complex.Sqrt(discriminant), 1.0f / 3.0f);

            roots[0] = u + v - b / (3 * a);
            roots[1] = -(u + v) / 2 - b / (3 * a) + Complex.Sqrt(3) * (u - v) / 2 * Complex.ImaginaryOne;
            roots[2] = -(u + v) / 2 - b / (3 * a) - Complex.Sqrt(3) * (u - v) / 2 * Complex.ImaginaryOne;
        }
        else
        {
            Complex A = Complex.Pow(-q / 2 + Complex.Sqrt(discriminant), 1.0f / 3.0f);
            Complex B = Complex.Pow(-q / 2 - Complex.Sqrt(discriminant), 1.0f / 3.0f);

            roots[0] = A + B - b / (3 * a);
            roots[1] = -(A + B) / 2 - b / (3 * a) + Complex.Sqrt(3) * (A - B) / 2 * omega;
            roots[2] = -(A + B) / 2 - b / (3 * a) - Complex.Sqrt(3) * (A - B) / 2 * omega * omega;
        }

        return roots;
    }

    // Function to check if a number is a real number (does not have an imaginary component)
    public static bool IsRealNumber(Complex number)
    {
        return Math.Abs(number.Imaginary) < double.Epsilon;
    }

    public static float setSmallNumberToZero(float number)
    {
        const float threshold = 0.00001f;

        if (Mathf.Abs(number) <= threshold)
        {
            return 0f;
        }
        else
        {
            return number;
        }
    }

    public static float[] calcCubicEquationsRealRoots(float a1, float b, float c, float d)  // solve cubic equation according to cardano
    {
        float p, q, u, v;
        float r, alpha;
        float[] roots = new float[3];
        if (a1 != 0f)
        {
            float a = b / a1;
            b = c / a1;
            c = d / a1;

            p = -(a * a / 3.0f) + b;
            q = (2.0f / 27.0f * a * a * a) - (a * b / 3.0f) + c;
            d = q * q / 4.0f + p * p * p / 27.0f;
            if (Math.Abs(d) < Math.Pow(10.0f, -11.0f))
                d = 0f;
            // 3 cases D > 0, D == 0 and D < 0
            if (d > 1e-20f)
            {
                u = Xroot(-q / 2.0f + Mathf.Sqrt(d), 3.0f);
                v = Xroot(-q / 2.0f - Mathf.Sqrt(d), 3.0f);
                roots[0] = u + v - a / 3.0f;
            }
            if (Math.Abs(d) <= 1e-20f)
            {
                u = Xroot(-q / 2.0f, 3.0f);
                v = Xroot(-q / 2.0f, 3.0f);
                roots[0] = u + v - a / 3.0f;
                roots[1] = -(u + v) / 2.0f - a / 3.0f;
            }
            if (d < -1e-20f)
            {
                r = Mathf.Sqrt(-p * p * p / 27.0f);
                alpha = Mathf.Atan(Mathf.Sqrt(-d) / q * 2.0f);
                if (q > 0f)                         // if q > 0 the angle becomes  PI + alpha
                    alpha = Mathf.PI + alpha;

                roots[0] = Xroot(r, 3.0f) * (Mathf.Cos((6.0f * Mathf.PI - alpha) / 3.0f) + Mathf.Cos(alpha / 3.0f)) - a / 3.0f;
                roots[1] = Xroot(r, 3.0f) * (Mathf.Cos((2.0f * Mathf.PI + alpha) / 3.0f) + Mathf.Cos((4.0f * Mathf.PI - alpha) / 3.0f)) - a / 3.0f;
                roots[2] = Xroot(r, 3.0f) * (Mathf.Cos((4.0f * Mathf.PI + alpha) / 3.0f) + Mathf.Cos((2.0f * Mathf.PI - alpha) / 3.0f)) - a / 3.0f;
            }
            return roots;
        }
        return roots;
    }

    public static float GetRoot (float r0, float r1, float r2, float z0, float z1, float z2, float g, float maxIterations)
    {
        float n0 = r0*z0;
        float n1 = r1*z1;
        float n2 = r2*z2;
        float t0 = z2 - 1f;
        if(r0 == 1f)
        {
            t0 = z0 - 1f;
        }
        if(r1 == 1f)
        {
            t0 = z1 - 1f;
        }
        if(r2 == 1f)
        {
            t0 = z2 - 1f;
        }
    
        float t1;
        if(g < 0)
        {
            t1 = 0;
        }
        else
        {
            t1 = -1 + new Vector3(n0, n1, n2).magnitude;
        }
        float t = 0;
        for(int i = 0; i<maxIterations; i++)
        {
            t = (t0 + t1)/2f;
            if(t == t0 || t == t1) {break;}
            float ratio0 = n0/(t+r0);
            float ratio1 = n1/(t+r1);
            float ratio2 = n2/(t+r2);
            g = ratio0*ratio0 + ratio1*ratio1 + ratio2*ratio2 - 1f;
            if(g>0) {t0 = t;}
            else
            {
                if(g<0) {t1 = t;}
                else {break;}
            }
        }

        return t;
    }

    public static Color[] InterpolateColor(float[] aInput, float min, float max)
    {
        Color[] colors = new Color[aInput.Length];
        for(int i = 0; i<aInput.Length; i++)
        {
            float range = max - min;
            float f = 2f * (aInput[i] - min) / range; // float in the range of 0 - 2
            if (f > 1)
                colors[i] = Color.Lerp(Color.yellow, Color.red, f - 1f); // top half
            colors[i] = Color.Lerp(Color.blue, Color.yellow, f); // bottom half
        }
        return colors;
    }

    public static void eigen(float[,] matrix , out Vector2[] vectors, out float[] values)
    {
        float b = matrix[1,1] + matrix[0,0];
        float c = matrix[1,1]*matrix[0,0] - matrix[0,1]*matrix[1,0];
        float val1 = 0.5f*(b + Mathf.Sqrt(b*b-4f*c));
        float val2 = 0.5f*(b - Mathf.Sqrt(b*b-4f*c));
        vectors = new Vector2[2];
        values = new float[2];


        if(val1>val2)
        {
            float val1_ = val1;
            val1 = val2;
            val2 = val1_; 
        }
        vectors[0] = new Vector2(1f, -matrix[0,1]/(matrix[0,0]-val1));
        vectors[1] = new Vector2(1f, -matrix[0,1]/(matrix[0,0]-val2));

        vectors[0] = vectors[0].normalized;
        vectors[1] = vectors[1].normalized;

        values[0] = val1;
        values[1] = val2;
    }

    public static float[,] matmul(float[,] m1, float[,] m2)
    {
        float[,] mr = new float[m1.GetLength(0), m2.GetLength(1)];
        for (int i = 0; i < m1.GetLength(0); i++)
        {
            for (int j = 0; j < m2.GetLength(1); j++)
            {
                for (int l = 0; l < m1.GetLength(1); l++)
                {
                    mr[i,j] += m1[i,l]*m2[l,j];
                }
            }
        }
        return mr;
    }

    public static float calDetMatrix3x3(float[,] matrix)
    {
        float det = 0; 
        for (int i = 0; i < 3; i++)
            det = det + (matrix[0, i] * (matrix[1, (i + 1) % 3] * matrix[2, (i + 2) % 3] - matrix[1, (i + 2) % 3] * matrix[2, (i + 1) % 3]));
        return det;
    }

    public static float calDetMatrix4x4(float[,] matrix)
    {
        float[,] matrix11 = new float[3, 3] { { matrix[1, 1], matrix[1, 2], matrix[1, 3] }, { matrix[2, 1], matrix[2, 2], matrix[2, 3] }, 
            { matrix[3, 1], matrix[3, 2], matrix[3, 3] } };
        float[,] matrix21 = new float[3, 3] { { matrix[0, 1], matrix[0, 2], matrix[0, 3] }, { matrix[2, 1], matrix[2, 2], matrix[2, 3] },
            { matrix[3, 1], matrix[3, 2], matrix[3, 3] } };
        float[,] matrix31 = new float[3, 3] { { matrix[0, 1], matrix[0, 2], matrix[0, 3] }, { matrix[1, 1], matrix[1, 2], matrix[1, 3] },
            { matrix[3, 1], matrix[3, 2], matrix[3, 3] } };
        float[,] matrix41 = new float[3, 3] { { matrix[0, 1], matrix[0, 2], matrix[0, 3] }, { matrix[1, 1], matrix[1, 2], matrix[1, 3] },
            { matrix[2, 1], matrix[2, 2], matrix[2, 3] } };
        return matrix[0, 0] * calDetMatrix3x3(matrix11) - matrix[1, 0] * calDetMatrix3x3(matrix21) + matrix[2, 0] * calDetMatrix3x3(matrix31)
            - matrix[3, 0] * calDetMatrix3x3(matrix41);
    }

    // Method to find the closest points on two lines A and B
    public static void FindClosestPoints(Vector3 A1, Vector3 A2, Vector3 B1, Vector3 B2, out Vector3 closestPointLineA, out Vector3 closestPointLineB)
    {
        Vector3 directionA = A2 - A1;
        Vector3 directionB = B2 - B1;
        Vector3 crossProduct = Vector3.Cross(directionA, directionB);

        // Check if lines are parallel
        if (crossProduct.sqrMagnitude < 0.001f)
        {
            closestPointLineA = A1;
            closestPointLineB = B1;
            return;
        }

        Vector3 lineVector = A1 - B1;
        float a = Vector3.Dot(directionA, directionA);
        float b = Vector3.Dot(directionA, directionB);
        float c = Vector3.Dot(directionB, directionB);
        float d = Vector3.Dot(directionA, lineVector);
        float e = Vector3.Dot(directionB, lineVector);

        float denominator = a * c - b * b;

        // Calculate the parameters for the lines
        float s = (b * e - c * d) / denominator;
        float t = (a * e - b * d) / denominator;

        closestPointLineA = A1 + directionA * s;
        closestPointLineB = B1 + directionB * t;        
    }

    public static float getRotationX(GameObject obj)
    {
        float rotation = obj.transform.eulerAngles.x;
        return rotation <= 180f ? rotation : rotation - 360f;
    }

    public static float getRotationY(GameObject obj)
    {
        float rotation = obj.transform.eulerAngles.y;
        return rotation <= 180f ? rotation : rotation - 360f;
    }

    public static float getRotationZ(GameObject obj)
    {
        float rotation = obj.transform.eulerAngles.z;
        return rotation <= 180f ? rotation : rotation - 360f;
    }

    public static void DrawPlane(Vector3 position, Vector3 normal)
    {
        Vector3 v3;

        if (normal.normalized != Vector3.forward)
            v3 = Vector3.Cross(normal, Vector3.forward).normalized * normal.magnitude;
        else
            v3 = Vector3.Cross(normal, Vector3.up).normalized * normal.magnitude;

        var corner0 = position + v3;
        var corner2 = position - v3;

        var q = Quaternion.AngleAxis(90f, normal);
        v3 = q * v3;
        var corner1 = position + v3;
        var corner3 = position - v3;

        Debug.DrawLine(corner0, corner2, Color.green);
        Debug.DrawLine(corner1, corner3, Color.green);
        Debug.DrawLine(corner0, corner1, Color.green);
        Debug.DrawLine(corner1, corner2, Color.green);
        Debug.DrawLine(corner2, corner3, Color.green);
        Debug.DrawLine(corner3, corner0, Color.green);
        Debug.DrawRay(position, normal, Color.red);
    }

    // Function to return whether the middle plane is between 2 other parallel planes
    public static bool IsPlaneBetween(Vector3 plane1Center, Vector3 middlePlaneCenter, Vector3 plane2Center)
    {
        // Define vectors along the common axis
        Vector3 axis = (plane2Center - plane1Center).normalized;
        Vector3 toMiddle = middlePlaneCenter - plane1Center;

        // Project the toMiddle vector onto the axis
        float projection = Vector3.Dot(toMiddle, axis);

        // Check if the middle plane is between the other two along the common axis
        return projection > 0 && projection < Vector3.Distance(plane1Center, plane2Center);
    }

    public static double[] SAT(Vector3 axis, Vector3[] corners, double minAlong, double maxAlong)
    {

        minAlong = double.MaxValue;
        maxAlong = -1 * double.MaxValue;

        for (int i = 0; i < corners.Length; i++)
        {

            double dotVal = Vector3.Dot(corners[i], axis);
            if (dotVal < minAlong) minAlong = dotVal;
            if (dotVal > maxAlong) maxAlong = dotVal;

        }
        double[] r = { minAlong, maxAlong };
        return r;

    }

    public static bool overlaps(double min1, double max1, double min2, double max2)
    {
        return isBetweenOrdered(min2, min1, max1) || isBetweenOrdered(min1, min2, max2);
    }

    public static bool isBetweenOrdered(double val, double lowerBound, double upperBound)
    {
        return lowerBound <= val && val <= upperBound;
    }

}
public static class Vector2Extensions
{
    public static Vector2 Rotate(this Vector2 vector, float angle)
    {
        // Create a Quaternion that represents the desired rotation
        Quaternion rotation = Quaternion.Euler(0, 0, angle);

        // Rotate the Vector3
        Vector3 rotatedV3 = rotation * vector;

        // Convert the rotated Vector3 back to Vector2
        return new Vector2(rotatedV3.x, rotatedV3.y);
    }
}

