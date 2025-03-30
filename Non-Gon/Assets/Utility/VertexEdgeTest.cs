using System;
using System.Collections.Generic;
using Complex = System.Numerics.Complex;
using UnityEngine;

public class VertexEdgeTest : MonoBehaviour
{
    public static GameObject cylinder1Global;
    public static GameObject cylinder2Global;
    public static float s1Global;
    public static float s2Global;
    public static float r1Global;
    public static float r2Global;
    public static float aGlobal;
    public static float bGlobal;
    public static float cGlobal;
    public static float alphaGlobal;
    public static Vector3 commonNormalGlobal;

    // Perform Vertex-Edge Test
    public static bool VertexEdgeTestFunction(float s1, float s2, float r1, float r2, float a, float b, float c, float alpha, GameObject cylinder1, GameObject cylinder2, Vector3 commonNormal)
    {
        Debug.Log("Vertex-Edge Test");
        cylinder1Global = cylinder1;
        cylinder2Global = cylinder2;
        s1Global = s1;
        s2Global = s2;
        r1Global = r1;
        r2Global = r2;
        aGlobal = a;
        bGlobal = b;
        cGlobal = c;
        alphaGlobal = alpha * Mathf.Deg2Rad;
        commonNormalGlobal = commonNormal;

        Vector2[] cylinder1Vertices = { new Vector2(1f,1f), new Vector2(-1f, 1f), new Vector2(-1f, -1f), new Vector2(1f, -1f)};
        Vector2[] cylinder2Vertices = { new Vector2(1f,1f), new Vector2(-1f, 1f), new Vector2(-1f, -1f), new Vector2(1f, -1f)};

        //Debug.Log("---------------------");
        int vertexCylinder1 = 0;
        int vertexCylinder2 = 0;
        foreach (Vector2 cylinder1Vertex in cylinder1Vertices) {
            foreach (Vector2 cylinder2Vertex in cylinder2Vertices)
            {
                Vector2 v1AndV2 = ApplyTransformation(cylinder1Vertex, cylinder2Vertex);
                float u1 = a + Mathf.Sqrt(Mathf.Abs(Mathf.Pow(r2, 2f) - Mathf.Pow(v1AndV2.y, 2f)));
                float u2 = a - Mathf.Sqrt(Mathf.Abs(Mathf.Pow(r2, 2f) - Mathf.Pow(v1AndV2.y, 2f)));
                float w1 = u1 - a;
                float w2 = u2 - a;
                float v1 = v1AndV2.x;
                float v2 = v1AndV2.y;
                //Debug.Log("u1: " + u1);
                //Debug.Log("u2: " + u2);
                //Debug.Log("v1: " + v1);
                //Debug.Log("v2: " + v2);
                if (!float.IsNaN(u1)) {
                    if ((u1 >= a - Mathf.Min(r1, r2)) & (u1 <= Mathf.Max(r1, r2)))
                    {
                        if (checkVertices(vertexCylinder1, vertexCylinder2, u1, w1, v1, v2))
                        {
                            return true;
                        }
                    }
                }
                if (!float.IsNaN(u2))
                {
                    if ((u2 >= a - Mathf.Min(r1, r2)) & (u2 <= Mathf.Max(r1, r2)))
                    {   
                        if (checkVertices(vertexCylinder1, vertexCylinder2, u2, w2, v1, v2))
                        {
                            return true;
                        }
                    }
                }
                vertexCylinder2++;
            }
            vertexCylinder1++;
            vertexCylinder2 = 0;
        }
        //Debug.Log("---------------------");
        // If no intersection found
        return false;
    }

    public static bool rectangleIntersection(float u, float w)
    {
        Debug.Log("Rectangle Intersection");
        Vector3 cylinder1Position = cylinder1Global.transform.position;
        Vector3 cylinder2Position = cylinder2Global.transform.position;

        Vector3 Circle1Center1 = cylinder1Position + cylinder1Global.transform.forward * s1Global;
        Vector3 Circle2Center1 = cylinder1Position - cylinder1Global.transform.forward * s1Global;
        Vector3 Circle1Center2 = cylinder2Position + cylinder2Global.transform.forward * s2Global;
        Vector3 Circle2Center2 = cylinder2Position - cylinder2Global.transform.forward * s2Global;
        Vector3 PlanePosition1 = cylinder1Position + commonNormalGlobal * u;
        Vector3 PlanePosition2 = cylinder2Position + commonNormalGlobal * w;

        //Debug.DrawLine(cylinder2Position, PlanePosition2, Color.red);

        Vector3[] point1 = CirclePlaneIntersection.FindIntersectionPoints(PlanePosition1, commonNormalGlobal, Circle1Center1, cylinder1Global.transform.forward, r1Global);
        Vector3[] point2 = CirclePlaneIntersection.FindIntersectionPoints(PlanePosition1, commonNormalGlobal, Circle2Center1, cylinder1Global.transform.forward, r1Global);
        Vector3[] point3 = CirclePlaneIntersection.FindIntersectionPoints(PlanePosition2, commonNormalGlobal, Circle1Center2, cylinder2Global.transform.forward, r2Global);
        Vector3[] point4 = CirclePlaneIntersection.FindIntersectionPoints(PlanePosition2, commonNormalGlobal, Circle2Center2, cylinder2Global.transform.forward, r2Global);

        if (point1.Length == 2 & point2.Length == 2 & point3.Length == 2 & point4.Length == 2)
        {
            List<Vector3> Q1Vertices = new List<Vector3>();
            if (Vector3.Distance(point1[1], point2[0]) < Vector3.Distance(point1[1], point2[1]))
            {
                Q1Vertices.Add(point1[0]);
                Q1Vertices.Add(point1[1]);
                Q1Vertices.Add(point2[0]);
                Q1Vertices.Add(point2[1]);
            }
            else
            {
                Q1Vertices.Add(point1[0]);
                Q1Vertices.Add(point1[1]);
                Q1Vertices.Add(point2[1]);
                Q1Vertices.Add(point2[0]);
            }
            List<Vector3> Q2Vertices = new List<Vector3>();
            if (Vector3.Distance(point3[1], point4[0]) < Vector3.Distance(point3[1], point4[1]))
            {
                Q2Vertices.Add(point3[0]);
                Q2Vertices.Add(point3[1]);
                Q2Vertices.Add(point4[0]);
                Q2Vertices.Add(point4[1]);
            }
            else
            {
                Q2Vertices.Add(point3[0]);
                Q2Vertices.Add(point3[1]);
                Q2Vertices.Add(point4[1]);
                Q2Vertices.Add(point4[0]);
            }
            
            /*
            Debug.DrawLine(Q1Vertices[0], Q1Vertices[1], Color.blue);
            Debug.DrawLine(Q1Vertices[1], Q1Vertices[2], Color.blue);
            Debug.DrawLine(Q1Vertices[2], Q1Vertices[3], Color.blue);
            Debug.DrawLine(Q1Vertices[3], Q1Vertices[0], Color.blue);

            Debug.DrawLine(Q2Vertices[0], Q2Vertices[1], Color.blue);
            Debug.DrawLine(Q2Vertices[1], Q2Vertices[2], Color.blue);
            Debug.DrawLine(Q2Vertices[2], Q2Vertices[3], Color.blue);
            Debug.DrawLine(Q2Vertices[3], Q2Vertices[0], Color.blue);
            */

            // Separating Axis Test (SAT) in Vertex Edge Test - Check for overlap along each axis of Q1 and Q2
            if (RectanglesIntersection.RectanglesIntersect(Q1Vertices, Q2Vertices))
            {
                Debug.Log("Intersection detected by SAT in Vertex Edge Test");
                return true;
            }
        }
        return false;
    }

    public static bool checkVertices(int vertexCylinder1, int vertexCylinder2, float u, float w, float v1, float v2)
    {
        Debug.Log("Check Vertices");
        // Case 1: top and bottom verification
        // Case 2 is in SideVerifications function

        // Vertex K1 Edge K2 L2
        if (vertexCylinder1 == 0 & vertexCylinder2 == 0)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex K1 Edge M2 N2
        if (vertexCylinder1 == 0 & vertexCylinder2 == 2)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex K1 Edge L2 M2
        if (vertexCylinder1 == 0 & vertexCylinder2 == 1)
        {
            float f = 1f;
            float g = Mathf.Cos(alphaGlobal);
            float h = (cGlobal + s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }
        // Vertex K1 Edge K2 N2
        if (vertexCylinder1 == 0 & vertexCylinder2 == 3)
        {
            float f = 1f;
            float g = -Mathf.Cos(alphaGlobal);
            float h = (cGlobal - s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }

        // Vertex L1 Edge K2 L2
        if (vertexCylinder1 == 1 & vertexCylinder2 == 0)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex L1 Edge M2 N2
        if (vertexCylinder1 == 1 & vertexCylinder2 == 2)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex L1 Edge L2 M2
        if (vertexCylinder1 == 1 & vertexCylinder2 == 1)
        {
            float f = -1f;
            float g = Mathf.Cos(alphaGlobal);
            float h = (cGlobal + s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }
        // Vertex L1 Edge K2 N2
        if (vertexCylinder1 == 1 & vertexCylinder2 == 3)
        {
            float f = -1f;
            float g = -Mathf.Cos(alphaGlobal);
            float h = (cGlobal - s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }

        // Vertex M1 Edge K2 L2
        if (vertexCylinder1 == 2 & vertexCylinder2 == 0)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex M1 Edge M2 N2
        if (vertexCylinder1 == 2 & vertexCylinder2 == 2)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex M1 Edge L2 M2
        if (vertexCylinder1 == 2 & vertexCylinder2 == 1)
        {
            float f = -1f;
            float g = Mathf.Cos(alphaGlobal);
            float h = (cGlobal + s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }
        // Vertex M1 Edge K2 N2
        if (vertexCylinder1 == 2 & vertexCylinder2 == 3)
        {
            float f = -1f;
            float g = -Mathf.Cos(alphaGlobal);
            float h = (cGlobal - s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }

        // Vertex N1 Edge K2 L2
        if (vertexCylinder1 == 3 & vertexCylinder2 == 0)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex N1 Edge M2 N2
        if (vertexCylinder1 == 3 & vertexCylinder2 == 2)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex N1 Edge L2 M2
        if (vertexCylinder1 == 3 & vertexCylinder2 == 1)
        {
            float f = 1f;
            float g = Mathf.Cos(alphaGlobal);
            float h = (cGlobal + s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }
        // Vertex N1 Edge K2 N2
        if (vertexCylinder1 == 3 & vertexCylinder2 == 3)
        {
            float f = 1f;
            float g = -Mathf.Cos(alphaGlobal);
            float h = (cGlobal - s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }

        // Vertex K2 Edge K1 L1
        if (vertexCylinder1 == 0 & vertexCylinder2 == 0)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex K2 Edge M1 N1
        if (vertexCylinder1 == 2 & vertexCylinder2 == 0)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex K2 Edge L1 M1
        if (vertexCylinder1 == 1 & vertexCylinder2 == 0)
        {
            float f = -1f;
            float g = -Mathf.Cos(alphaGlobal);
            float h = (cGlobal + s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }
        // Vertex K2 Edge K1 N1
        if (vertexCylinder1 == 3 & vertexCylinder2 == 0)
        {
            float f = 1f;
            float g = -Mathf.Cos(alphaGlobal);
            float h = (cGlobal + s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }

        // Vertex L2 Edge K1 L1
        if (vertexCylinder1 == 0 & vertexCylinder2 == 1)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex L2 Edge M1 N1
        if (vertexCylinder1 == 2 & vertexCylinder2 == 1)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex L2 Edge L1 M1
        if (vertexCylinder1 == 1 & vertexCylinder2 == 1)
        {
            float f = -1f;
            float g = Mathf.Cos(alphaGlobal);
            float h = (cGlobal + s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }
        // Vertex L2 Edge K1 N1
        if (vertexCylinder1 == 3 & vertexCylinder2 == 1)
        {
            float f = 1f;
            float g = Mathf.Cos(alphaGlobal);
            float h = (cGlobal + s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }

        // Vertex M2 Edge K1 L1
        if (vertexCylinder1 == 0 & vertexCylinder2 == 2)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex M2 Edge M1 N1
        if (vertexCylinder1 == 2 & vertexCylinder2 == 2)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex M2 Edge L1 M1
        if (vertexCylinder1 == 1 & vertexCylinder2 == 2)
        {
            float f = -1f;
            float g = Mathf.Cos(alphaGlobal);
            float h = (cGlobal - s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }
        // Vertex M2 Edge K1 N1
        if (vertexCylinder1 == 3 & vertexCylinder2 == 2)
        {
            float f = 1f;
            float g = Mathf.Cos(alphaGlobal);
            float h = (cGlobal - s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }

        // Vertex N2 Edge K1 L1
        if (vertexCylinder1 == 0 & vertexCylinder2 == 3)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex N2 Edge M1 N1
        if (vertexCylinder1 == 2 & vertexCylinder2 == 3)
        {
            if (rectangleIntersection(u, w))
            {
                return true;
            }
        }
        // Vertex N2 Edge L1 M1
        if (vertexCylinder1 == 1 & vertexCylinder2 == 3)
        {
            float f = -1f;
            float g = -Mathf.Cos(alphaGlobal);
            float h = (cGlobal - s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }
        // Vertex N2 Edge K1 N1
        if (vertexCylinder1 == 3 & vertexCylinder2 == 3)
        {
            float f = 1f;
            float g = -Mathf.Cos(alphaGlobal);
            float h = (cGlobal - s2Global) * Mathf.Sin(alphaGlobal);
            if (SideVerification(f, g, h, u, w))
            {
                return true;
            }
        }
        return false;
    }

    // Case 2: side verification
    public static bool SideVerification(float f, float g, float h, float u, float w)
    {
        Debug.Log("Side Verification");

        // Quartic equation coefficients
        float lambda1 = -Mathf.Pow(f, 4f) + 2f * Mathf.Pow(f, 2f) * Mathf.Pow(g, 2f) - Mathf.Pow(g, 4f);
        float lambda2 = -4f * aGlobal * Mathf.Pow(f, 2f) * Mathf.Pow(g, 2f) + 4f * aGlobal * Mathf.Pow(g, 4f);
        float lambda3 = 2f * Mathf.Pow(r1Global, 2f) * Mathf.Pow(f, 4f) + 2f * Mathf.Pow(aGlobal, 2f) * Mathf.Pow(f, 2f) * Mathf.Pow(g, 2f) - 2f * Mathf.Pow(r1Global, 2f) * Mathf.Pow(f, 2f) * Mathf.Pow(g, 2f) -
            2f * Mathf.Pow(r1Global, 2f) * Mathf.Pow(f, 2f) * Mathf.Pow(g, 2f) - 6f * Mathf.Pow(aGlobal, 2f) * Mathf.Pow(g, 4f) + 2f * Mathf.Pow(r2Global, 2f) * Mathf.Pow(g, 4f) -
            2f * Mathf.Pow(f, 2f) * Mathf.Pow(h, 2f) - 2f * Mathf.Pow(g, 2f) * Mathf.Pow(h, 2f);
        float lambda4 = 4f * aGlobal * Mathf.Pow(r1Global, 2f) * Mathf.Pow(f, 2f) * Mathf.Pow(g, 2f) + 4f * Mathf.Pow(aGlobal, 3f) * Mathf.Pow(g, 4f) - 4f * aGlobal * Mathf.Pow(r2Global, 2f) * Mathf.Pow(g, 4f) +
            4f * aGlobal * Mathf.Pow(g, 2f) * Mathf.Pow(h, 2f);
        float lambda5 = -Mathf.Pow(r1Global, 4f) * Mathf.Pow(f, 4f) - 2f * Mathf.Pow(aGlobal, 2f) * Mathf.Pow(r1Global, 2f) * Mathf.Pow(f, 2f) * Mathf.Pow(g, 2f) + 2f * Mathf.Pow(r1Global, 2f) * Mathf.Pow(r2Global, 2f) *
            Mathf.Pow(f, 2f) * Mathf.Pow(g, 2f) - Mathf.Pow(aGlobal, 4f) * Mathf.Pow(g, 4f) + 2f * Mathf.Pow(aGlobal, 2f) * Mathf.Pow(r2Global, 2f) * Mathf.Pow(g, 4f) - Mathf.Pow(r2Global, 4f) * Mathf.Pow(g, 4f) +
            2f * Mathf.Pow(r1Global, 2f) * Mathf.Pow(f, 2f) * Mathf.Pow(h, 2f) - 2f * Mathf.Pow(aGlobal, 2f) * Mathf.Pow(g, 2f) * Mathf.Pow(h, 2f) + 2f * Mathf.Pow(r2Global, 2f) * Mathf.Pow(g, 2f) * Mathf.Pow(h, 2f) -
            Mathf.Pow(h, 4f);

        // Solve quartic equation for u
        List<float> solutions = SolveRealQuartic(lambda1, lambda2, lambda3, lambda4, lambda5);

        if (solutions.Count == 2f)
        {
            if (((solutions[0] >= aGlobal - Mathf.Min(r1Global, r2Global)) & (solutions[0] <= aGlobal + Mathf.Max(r1Global, r2Global))) || 
                ((solutions[1] >= aGlobal - Mathf.Min(r1Global, r2Global)) & (solutions[1] <= aGlobal + Mathf.Max(r1Global, r2Global))))
            {
                if (rectangleIntersection(u, w))
                {
                    return true;
                }
            }
        }
        else if (solutions.Count == 1f)
        {
            if ((solutions[0] >= aGlobal - Mathf.Min(r1Global, r2Global)) & (solutions[0] <= aGlobal + Mathf.Max(r1Global, r2Global)))
            {
                if (rectangleIntersection(u, w))
                {
                    return true;
                }
            }
        }
        return false;
    }

    // Function to get the rectangle edges from a set of points
   

    public static List<float> SolveRealQuartic(float a, float b, float c, float d, float e)
    {
        Debug.Log("Solve Real Quartic");
        double D0 = c * c - 3 * b * d + 12 * a * e;
        double D1 = 2 * c * c * c - 9 * b * c * d + 27 * b * b * e + 27 * a * d * d - 72 * a * c * e;
        double p = (8 * a * c - 3 * b * b) / (8 * a * a);
        double q = (b * b * b - 4 * a * b * c + 8 * a * a * d) / (8 * a * a * a);
        Complex Q = Complex.Pow((D1 + Complex.Sqrt(D1 * D1 - 4 * D0 * D0 * D0)) / 2, 1.0 / 3.0);
        Complex S = Complex.Sqrt(-2 * p / 3 + (Q + D0 / Q) / (3 * a)) / 2;
        Complex u = Complex.Sqrt(-4 * S * S - 2 * p + q / S) / 2;
        Complex v = Complex.Sqrt(-4 * S * S - 2 * p - q / S) / 2;
        Complex x1 = -b / (4 * a) - S + u;
        Complex x2 = -b / (4 * a) - S - u;
        Complex x3 = -b / (4 * a) + S + v;
        Complex x4 = -b / (4 * a) + S - v;
        List<float> roots = new List<float>();
        int count = 0;
        if (x1.Imaginary == 0f)
        {
            roots.Add((float)x1.Real);
            count++;
        }
        if (x2.Imaginary == 0f)
        {
            roots.Add((float)x2.Real);
            count++;
        }
        if (x3.Imaginary == 0f)
        {
            roots.Add((float)x3.Real);
            count++;
        }
        if (x4.Imaginary == 0f)
        {
            roots.Add((float)x4.Real);
            count++;
        }
        //Debug.Log(count);
        return roots;
    }

    // Function to apply the transformation [T_B]_A to a vertex
    private static Vector2 ApplyTransformation(Vector2 v1Ands1Signs, Vector2 v2Ands2Signs)
    {
        Debug.Log("Apply Transformation");
        float actualS1 = v1Ands1Signs.y * s1Global;
        float actualS2 = v2Ands2Signs.y * s2Global;
        /*
        Debug.Log(actualS1);
        Debug.Log(actualS2);
        Debug.Log(bGlobal);
        Debug.Log(cGlobal);
        Debug.Log(alphaGlobal);
        Debug.Log("----------------");
        */
        // (s1 - b - (c - s2) cos alpha) / sin alpha
        float v2 = v2Ands2Signs.x * ((actualS1 - bGlobal - cGlobal * Mathf.Cos(alphaGlobal) - actualS2 * Mathf.Cos(alphaGlobal)) / Mathf.Sin(alphaGlobal));
        // (s1 - b - (c - s2) cos alpha) / sin alpha
        float v1 = v1Ands1Signs.x * (v2 * Mathf.Cos(alphaGlobal) - actualS2 * Mathf.Sin(alphaGlobal) - cGlobal * Mathf.Sin(alphaGlobal));

        return new Vector2(v1, v2);
    }

}