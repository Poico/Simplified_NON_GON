using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using static Utility;
using static ProximityQuery2D;
using Complex = System.Numerics.Complex;

public class ProximityQuery3D : MonoBehaviour
{
    public static bool pass = false;
    //public static bool PointSphere3D(Vector3 positionPoint, Vector3 positionSphere, double sphereRadius)
    //{
    //    Debug.Log("PointSphere3D");
    //    double pX = positionPoint.x;
    //    double pY = positionPoint.y;
    //    double pZ = positionPoint.z;
    //    double sX = positionSphere.x;
    //    double sY = positionSphere.y;
    //    double sZ = positionSphere.z;
    //
    //    double distance = Mathf.Sqrt((float)((pX - sX) * (pX - sX) +
    //                             (pY - sY) * (pY - sY) +
    //                             (pZ - sZ) * (pZ - sZ)));
    //    return distance < sphereRadius;
    //}

    //public static bool Point_AABB3D(Vector3 positionPoint, Vector3 positionAABB, double lengthAABB, double heightAABB, double widthAABB)
    //{
    //    Debug.Log("Point_AABB3D");
    //    double pX = positionPoint.x;
    //    double pY = positionPoint.y;
    //    double pZ = positionPoint.z;
    //    double minX = positionAABB.x - lengthAABB / 2;
    //    double maxX = positionAABB.x + lengthAABB / 2;
    //    double minY = positionAABB.y - widthAABB / 2;
    //    double maxY = positionAABB.y + widthAABB / 2;
    //    double minZ = positionAABB.z - heightAABB / 2;
    //    double maxZ = positionAABB.z + heightAABB / 2;
    //
    //    return (pX >= minX && pX <= maxX) &&
    //     (pY >= minY && pY <= maxY) &&
    //     (pZ >= minZ && pZ <= maxZ);
    //}

    //public static bool Sphere_AABB3D(Vector3 positionAABB, double lengthAABB, double heightAABB, double widthAABB, Vector3 positionSphere, double sphereRadius)
    //{
    //    Debug.Log("Sphere_AABB3D");
    //    double minX = positionAABB.x - lengthAABB / 2;
    //    double maxX = positionAABB.x + lengthAABB / 2;
    //    double minY = positionAABB.y - widthAABB / 2;
    //    double maxY = positionAABB.y + widthAABB / 2;
    //    double minZ = positionAABB.z - heightAABB / 2;
    //    double maxZ = positionAABB.z + heightAABB / 2;
    //    double sX = positionSphere.x;
    //    double sY = positionSphere.y;
    //    double sZ = positionSphere.z;
    //    double x = Mathf.Max((float)minX, Mathf.Min((float)sX, (float)maxX));
    //    double y = Mathf.Max((float)minY, Mathf.Min((float)sY, (float)maxY));
    //    double z = Mathf.Max((float)minZ, Mathf.Min((float)sZ, (float)maxZ));
    //
    //    double distance = Mathf.Sqrt((float)((x - sX) * (x - sX) + (y - sY) * (y - sY) + (z - sZ) * (z - sZ)));
    //
    //    return distance < sphereRadius;
    //}

    //public static bool SphereSphere3D(Vector3 positionSphere1, Vector3 positionSphere2, double sphereRadius1, double sphereRadius2)
    //{
    //    Debug.Log("SphereSphere3D");
    //    double sX1 = positionSphere1.x;
    //    double sX2 = positionSphere2.x;
    //    double sY1 = positionSphere1.y;
    //    double sY2 = positionSphere2.y;
    //    double sZ1 = positionSphere1.z;
    //    double sZ2 = positionSphere2.z;
    //
    //    double distance = Mathf.Sqrt((float)((sX1 - sX2) * (sX1 - sX2) +
    //                             (sY1 - sY2) * (sY1 - sY2) +
    //                             (sZ1 - sZ2) * (sZ1 - sZ2)));
    //    return distance < (sphereRadius1 + sphereRadius2);
    //}

    //public static bool AABB_AABB3D(Vector3 positionAABB1, Vector3 positionAABB2, double length1, double length2, double width1, double width2, double heigth1, double heigth2)
    //{
    //    Debug.Log("AABB_AABB3D");
    //    bool res = false;
    //    double x1 = positionAABB1.x;
    //    double y1 = positionAABB1.y;
    //    double x2 = positionAABB2.x;
    //    double y2 = positionAABB2.y;
    //    double z1 = positionAABB1.z;
    //    double z2 = positionAABB2.z;
    //    if ((x1 < (x2 + length2) && (x1 + length1) > x2) && (y1 < (y2 + width2) && (y1 + width1) > y2) && (z1 < (z2 + heigth2)) && ((z1 + heigth1) > z2))
    //    {
    //        res = true;
    //
    //    }
    //    return res;
    //}

    //public static bool OBB_OBB3D(Vector3[] normals, Vector3[] corners1, Vector3[] corners2)
    //{
    //    Debug.Log("OBB_OBB3D");
    //    for (int i = 0; i < normals.Length; i++)
    //    {
    //        float shape1Min = 0, shape1Max = 0, shape2Min = 0, shape2Max = 0;
    //        double[] shape1 = SAT(normals[i], corners1, shape1Min, shape1Max);
    //        double[] shape2 = SAT(normals[i], corners2, shape2Min, shape2Max);
    //        if (!Utility.overlaps(shape1[0], shape1[1], shape2[0], shape2[1]))
    //        {
    //            return false;
    //        }
    //    }
    //    return true;
    //}

    //For different sized boxes


    public static bool Cylinder_Cylinder_Chittawadigi(GameObject Cylinder1, GameObject Cylinder2)
    {
        
        GeometryCreator variables1 = Cylinder1.GetComponent<GeometryCreator>();
        GeometryCreator variables2 = Cylinder2.GetComponent<GeometryCreator>();
        float xradius1 = variables1.xradius;
        float xradius2 = variables2.xradius;
        float yradius1 = variables1.yradius;
        float yradius2 = variables2.yradius;
        // Get the positions and orientations of the cylinders
        Vector3 cylinder1Position = Cylinder1.transform.position;
        Vector3 cylinder2Position = Cylinder2.transform.position;

        // Get the Z-axis of each cylinder
        Vector3 zAxisCylinder1 = Cylinder1.transform.forward;
        Vector3 zAxisCylinder2 = Cylinder2.transform.forward;
        
        Vector3 closestPointAxisCylinder1;
        Vector3 closestPointAxisCylinder2;

        FindClosestPoints(cylinder1Position, cylinder1Position + zAxisCylinder1, cylinder2Position, cylinder2Position + zAxisCylinder2, out closestPointAxisCylinder1, out closestPointAxisCylinder2);

        // Find the common normal vector
        Vector3 commonNormal = closestPointAxisCylinder2 - closestPointAxisCylinder1;

        float b;
        float theta;
        float a;
        float alpha;
        float c;
        if (closestPointAxisCylinder1 == cylinder1Position & closestPointAxisCylinder2 == cylinder2Position)
        {
            b = 0f;
            //theta = Vector3.Angle(Cylinder1.transform.up, cylinder2Position - cylinder1Position);
            c = Mathf.Abs(Mathf.Round((cylinder2Position.z - cylinder1Position.z) * 100f)/100f);
            Vector3 endOfSecondRay = new Vector3(cylinder1Position.x, cylinder1Position.y, cylinder1Position.z + cylinder2Position.z - cylinder1Position.z);
            a = Mathf.Abs(Mathf.Round(Vector3.Distance(cylinder2Position, endOfSecondRay) * 100f)/100f);
            alpha = 0f;
            
        } else
        {
            // Calculating b
            b = Vector3.Distance(closestPointAxisCylinder1, cylinder1Position);

            //Debug.DrawRay(cylinder1Position, closestPointAxisCylinder1 - cylinder1Position, Color.yellow);

            // Calculating the end point of the first ray
            Vector3 endOfFirstRay = closestPointAxisCylinder1;

            //theta = Vector3.Angle(Cylinder1.transform.up, commonNormal);

            // Calculating a and alpha
            a = commonNormal.magnitude;

            //Debug.DrawRay(endOfFirstRay, commonNormal, Color.yellow);

            // Calculating the end point of the second ray
            Vector3 endOfSecondRay = endOfFirstRay + commonNormal;

            alpha = Mathf.Round(Vector3.Angle(Cylinder1.transform.forward, Cylinder2.transform.forward)*100f)/100f;

            if (alpha > 90f)
            {
                alpha = 180f - alpha;
            }
            // Calculating c
            c = Vector3.Distance(endOfSecondRay, cylinder2Position);

            //Debug.DrawRay(endOfSecondRay, cylinder2Position - endOfSecondRay, Color.yellow);
        }

        if (a < 0.00001f)
        {
            a = 0f;
        }
        if (b < 0.00001f)
        {
            b = 0f;
        }
        if (c < 0.00001f)
        {
            c = 0f;
        }
        if ((alpha < 0.00001f & alpha > 0f) || (alpha > -0.00001f & alpha < 0f))
        {
            alpha = 0f;
        }

        /*
        // Display the calculated values
        Debug.Log("b: " + b);
        Debug.Log("a: " + a);
        Debug.Log("Alpha: " + alpha);
        Debug.Log("c: " + c);
        */

        float s1 = yradius1 / 2f;
        float s2 = yradius2 / 2f;
        float r1 = xradius1;
        float r2 = xradius2;

        //Parallel Test
        if ((alpha == 0f || alpha == 180f) & b == 0f)
        {

            if (s1 + s2 >= Mathf.Abs(c) & r1 + r2 >= Mathf.Abs(a))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            //Non Parallel Test |b| <= s1 and |c| <= s2
            if (Mathf.Abs(b) <= s1 & Mathf.Abs(c) <= s2 & Mathf.Abs(a) <= r1 + r2)
            {
                return true;
            }
            else
            {
                bool doZAxisIntersect = false;
                // Check if lines intersect (closest points are identical)
                if (Vector3.Distance(closestPointAxisCylinder1, closestPointAxisCylinder2) <= 0.001f)
                {
                    Vector3 A1 = cylinder1Position;
                    Vector3 A2 = cylinder1Position + zAxisCylinder1;
                    Vector3 B1 = cylinder2Position;
                    Vector3 B2 = cylinder2Position + zAxisCylinder2;
                    Vector3 directionA = A2 - A1;
                    Vector3 directionB = B2 - B1;
                    Vector3 crossProduct = Vector3.Cross(directionA, directionB);
                    closestPointAxisCylinder2 = closestPointAxisCylinder2 + crossProduct;
                    commonNormal = closestPointAxisCylinder2 - closestPointAxisCylinder1;
                    doZAxisIntersect = true;
                }
                //DrawPlane(cylinder2Position, commonNormal);
                //DrawPlane(cylinder1Position, commonNormal);

                Vector3 Circle1Center1 = cylinder1Position + Cylinder1.transform.forward * s1;
                Vector3 Circle2Center1 = cylinder1Position - Cylinder1.transform.forward * s1;
                Vector3 Circle1Center2 = cylinder2Position + Cylinder2.transform.forward * s2;
                Vector3 Circle2Center2 = cylinder2Position - Cylinder2.transform.forward * s2;

                Vector3[] point1 = CirclePlaneIntersection.FindIntersectionPoints(cylinder1Position, commonNormal, Circle1Center1, Cylinder1.transform.forward, r1);
                Vector3[] point2 = CirclePlaneIntersection.FindIntersectionPoints(cylinder1Position, commonNormal, Circle2Center1, Cylinder1.transform.forward, r1);
                Vector3[] point3 = CirclePlaneIntersection.FindIntersectionPoints(cylinder2Position, commonNormal, Circle1Center2, Cylinder2.transform.forward, r2);
                Vector3[] point4 = CirclePlaneIntersection.FindIntersectionPoints(cylinder2Position, commonNormal, Circle2Center2, Cylinder2.transform.forward, r2);

                if (point1 == null || point2 == null || point3 == null || point4 == null)
                {
                    return false;
                }
                List<Vector3> Q1Vertices = new List<Vector3>();
                if (Vector3.Distance(point1[1], point2[0]) < Vector3.Distance(point1[1], point2[1]))
                {
                    Q1Vertices.Add(point1[0]);
                    Q1Vertices.Add(point1[1]);
                    Q1Vertices.Add(point2[0]);
                    Q1Vertices.Add(point2[1]);
                } else
                {
                    Q1Vertices.Add(point1[0]);
                    Q1Vertices.Add(point1[1]);
                    Q1Vertices.Add(point2[1]);
                    Q1Vertices.Add(point2[0]);
                }
                List<Vector3> Q2Vertices = new List<Vector3>();
                if (Vector3.Distance(point3[1], point4[0]) < Vector3.Distance(point3[1], point4[1]))
                {
                    if (doZAxisIntersect)
                    {
                        Q2Vertices.Add(point3[0]);
                        Q2Vertices.Add(point3[1]);
                        Q2Vertices.Add(point4[0]);
                        Q2Vertices.Add(point4[1]);
                    }
                    else
                    {
                        Q2Vertices.Add(point3[0] - commonNormal);
                        Q2Vertices.Add(point3[1] - commonNormal);
                        Q2Vertices.Add(point4[0] - commonNormal);
                        Q2Vertices.Add(point4[1] - commonNormal);
                    }
                }
                else
                {
                    if (doZAxisIntersect)
                    {
                        Q2Vertices.Add(point3[0]);
                        Q2Vertices.Add(point3[1]);
                        Q2Vertices.Add(point4[1]);
                        Q2Vertices.Add(point4[0]);
                    }
                    else
                    {
                        Q2Vertices.Add(point3[0] - commonNormal);
                        Q2Vertices.Add(point3[1] - commonNormal);
                        Q2Vertices.Add(point4[1] - commonNormal);
                        Q2Vertices.Add(point4[0] - commonNormal);
                    }
                }
                
                /*
                Debug.DrawLine(Q1Vertices[0], Q1Vertices[1], Color.red);
                Debug.DrawLine(Q1Vertices[1], Q1Vertices[2], Color.red);
                Debug.DrawLine(Q1Vertices[2], Q1Vertices[3], Color.red);
                Debug.DrawLine(Q1Vertices[3], Q1Vertices[0], Color.red);

                Debug.DrawLine(Q2Vertices[0], Q2Vertices[1], Color.red);
                Debug.DrawLine(Q2Vertices[1], Q2Vertices[2], Color.red);
                Debug.DrawLine(Q2Vertices[2], Q2Vertices[3], Color.red);
                Debug.DrawLine(Q2Vertices[3], Q2Vertices[0], Color.red);
                */
                // Separating Axis Test (SAT) - Check for overlap along each axis of Q1 and Q2
                if (!RectanglesIntersection.RectanglesIntersect(Q1Vertices, Q2Vertices))
                {
                    //Debug.Log("No intersection detected by SAT");
                    // No intersection detected by SAT
                    return false;
                }
                
                // Perform Vertex-Edge Test
                if (!VertexEdgeTest.VertexEdgeTestFunction(s1, s2, r1, r2, a, b, c, alpha, Cylinder1, Cylinder2, commonNormal))
                {
                    //Debug.Log("No intersection detected by Vertex-Edge Test");
                    // No intersection detected by Vertex-Edge Test
                    return false;
                } else {
                    // Intersection detected by Vertex-Edge Test
                    return true;
                }
            }
        }
        
    }
    
    public static float[] characteristicPolynomialEllipsoid(GameObject Ellipsoid1, GameObject Ellipsoid2)
    {
        GeometryCreator variables1 = Ellipsoid1.GetComponent<GeometryCreator>();
        float xradius1 = variables1.xradius; // Ellipsoid semi axis
        float yradius1 = variables1.yradius;
        float zradius1 = variables1.zradius;
        float xEllipsoid1 = Ellipsoid1.transform.position.x; // Ellipsoid position in space
        float yEllipsoid1 = Ellipsoid1.transform.position.y;
        float zEllipsoid1 = Ellipsoid1.transform.position.z;
        float alpha1 = Utility.getRotationX(Ellipsoid1) * Mathf.Deg2Rad;
        float beta1 = Utility.getRotationY(Ellipsoid1) * Mathf.Deg2Rad;
        float phi1 = Utility.getRotationZ(Ellipsoid1) * Mathf.Deg2Rad;
        float sinAlpha1 = Mathf.Sin(alpha1);
        float sinBeta1 = Mathf.Sin(beta1);
        float sinPhi1 = Mathf.Sin(phi1);
        float cosAlpha1 = Mathf.Cos(alpha1);
        float cosBeta1 = Mathf.Cos(beta1);
        float cosPhi1 = Mathf.Cos(phi1);

        GeometryCreator variables2 = Ellipsoid2.GetComponent<GeometryCreator>();
        float xradius2 = variables2.xradius; // Ellipsoid semi axis
        float yradius2 = variables2.yradius;
        float zradius2 = variables2.zradius;
        float xEllipsoid2 = Ellipsoid2.transform.position.x; // Ellipsoid position in space
        float yEllipsoid2 = Ellipsoid2.transform.position.y;
        float zEllipsoid2 = Ellipsoid2.transform.position.z;
        float alpha2 = Utility.getRotationX(Ellipsoid2) * Mathf.Deg2Rad;
        float beta2 = Utility.getRotationY(Ellipsoid2) * Mathf.Deg2Rad;
        float phi2 = Utility.getRotationZ(Ellipsoid2) * Mathf.Deg2Rad;
        float sinAlpha2 = Mathf.Sin(alpha2);
        float sinBeta2 = Mathf.Sin(beta2);
        float sinPhi2 = Mathf.Sin(phi2);
        float cosAlpha2 = Mathf.Cos(alpha2);
        float cosBeta2 = Mathf.Cos(beta2);
        float cosPhi2 = Mathf.Cos(phi2);

        float[,] matrixA = new float[4, 4];
        float[,] matrixB = new float[4, 4];

        float aa1 = xradius1 * xradius1;
        float bb1 = yradius1 * yradius1;
        float cc1 = zradius1 * zradius1;

        float aa2 = xradius2 * xradius2;
        float bb2 = yradius2 * yradius2;
        float cc2 = zradius2 * zradius2;

        //Quadric 1
        //x^2 (/aa1 /bb1 /cc1)
        float A1 = (cosBeta1 * cosBeta1 * cosPhi1 * cosPhi1) * bb1 * cc1 + 
            (cosBeta1 * cosBeta1 * sinPhi1 * sinPhi1) * aa1 * cc1 +
            (sinBeta1 * sinBeta1) * aa1 * bb1;
        //y^2
        float B1 = Mathf.Pow((sinBeta1 * sinAlpha1 * cosPhi1 + sinPhi1 * cosAlpha1), 2f) * bb1 * cc1 + 
            Mathf.Pow(sinPhi1 * sinBeta1 * sinAlpha1 + cosAlpha1 * cosPhi1, 2f) * aa1 * cc1 + 
            Mathf.Pow((sinAlpha1 * cosBeta1), 2f) * aa1 * bb1;
        //z^2
        float C1 = ((-sinBeta1 * cosAlpha1 * cosPhi1 + sinAlpha1 * sinPhi1) * (-sinBeta1 * sinAlpha1 * cosPhi1 + 
            sinAlpha1 * sinPhi1)) * bb1 * cc1 + 
            ((sinBeta1 * cosAlpha1 * sinPhi1 + cosPhi1 * sinAlpha1) * (sinBeta1 * cosAlpha1 * sinPhi1 + 
            cosPhi1 * sinAlpha1)) * aa1 * cc1 +
            (cosBeta1 * cosBeta1 * cosAlpha1 * cosAlpha1) * aa1 * bb1;
        //xy (must divide by 2f)
        float D1 = ((sinBeta1 * sinAlpha1 * cosPhi1 * cosBeta1 * cosPhi1 + 
            sinPhi1 * cosAlpha1 * cosBeta1 * cosPhi1) * bb1 * cc1 +
            (sinPhi1 * sinBeta1 * sinAlpha1 * cosBeta1 * sinPhi1 - cosAlpha1 * cosPhi1 * cosBeta1 * sinPhi1) * aa1 * cc1 +
            (-sinAlpha1 * cosBeta1 * sinBeta1) * aa1 * bb1);
        //yz (must divide by 2f)
        float E1 = (((-sinBeta1 * cosAlpha1 * cosPhi1 + sinAlpha1 * sinPhi1) * (sinBeta1 * sinAlpha1 * cosPhi1 + 
            sinPhi1 * cosAlpha1)) * bb1 * cc1 +
            ((sinBeta1 * cosAlpha1 * sinPhi1 + cosPhi1 * sinAlpha1) * (-sinPhi1 * sinBeta1 * sinAlpha1 + 
            cosAlpha1 * cosPhi1)) * aa1 * cc1 +
            (-cosBeta1 * cosAlpha1 * cosBeta1 * sinAlpha1) * aa1 * bb1);
        //zx (must divide by 2f)
        float F1 = ((-sinBeta1 * cosAlpha1 * cosPhi1 * cosBeta1 * cosPhi1 + sinAlpha1 * sinPhi1 * cosBeta1 * cosPhi1) * bb1 * cc1 +
            (-sinBeta1 + cosAlpha1 * sinPhi1 * cosBeta1 * sinPhi1 - cosPhi1 * sinAlpha1 * cosBeta1 * sinPhi1) * aa1 * cc1 +
            (cosBeta1 * cosAlpha1 * sinBeta1) * aa1 * bb1);

        //x (must divide by 2f)
        float G1 = - xEllipsoid1 * A1 - yEllipsoid1 * D1 - zEllipsoid1 * F1;
        //y (must divide by 2f)
        float H1 = - xEllipsoid1 * D1 - yEllipsoid1 * B1 - zEllipsoid1 * E1;
        //z (must divide by 2f)
        float I1 = - xEllipsoid1 * F1 - yEllipsoid1 * E1 - zEllipsoid1 * C1;
        //independent
        float J1 = (xEllipsoid1 * xEllipsoid1 * A1) + (xEllipsoid1 * yEllipsoid1 * 2f * D1) + (xEllipsoid1 * zEllipsoid1 * 2f * F1) + (yEllipsoid1 * yEllipsoid1 * B1) +
            (yEllipsoid1 * zEllipsoid1 * 2f * E1) + (zEllipsoid1 * zEllipsoid1 * C1) - (aa1 * bb1 * cc1);

        //Debug.Log(A1 + " " + B1 + " " + C1 + " " + D1 + " " + E1 + " " + F1 + " " + G1 + " " + H1 + " " + I1 + " " + J1);

        matrixA[0, 0] = A1;
        matrixA[0, 1] = D1;
        matrixA[0, 2] = F1;
        matrixA[0, 3] = G1;
        matrixA[1, 0] = D1;
        matrixA[1, 1] = B1;
        matrixA[1, 2] = E1;
        matrixA[1, 3] = H1;
        matrixA[2, 0] = F1;
        matrixA[2, 1] = E1;
        matrixA[2, 2] = C1;
        matrixA[2, 3] = I1;
        matrixA[3, 0] = G1;
        matrixA[3, 1] = H1;
        matrixA[3, 2] = I1;
        matrixA[3, 3] = J1;

        //Quadric 2
        //x^2
        float A2 = (cosBeta2 * cosBeta2 * cosPhi2 * cosPhi2) * bb2 * cc2 + 
            (cosBeta2 * cosBeta2 * sinPhi2 * sinPhi2) * aa2 * cc2 +
            (sinBeta2 * sinBeta2) * aa2 * bb2;
        //y^2
        float B2 = Mathf.Pow((sinBeta2 * sinAlpha2 * cosPhi2 + sinPhi2 * cosAlpha2), 2f) * bb2 * cc2 + 
            Mathf.Pow(sinPhi2 * sinBeta2 * sinAlpha2 + cosAlpha2 * cosPhi2, 2f) * aa2 * cc2 +
            Mathf.Pow((sinAlpha2 * cosBeta2), 2f) * aa2 * bb2;
        //z^2
        float C2 = ((-sinBeta2 * cosAlpha2 * cosPhi2 + sinAlpha2 * sinPhi2) * 
            (-sinBeta2 * sinAlpha2 * cosPhi2 + sinAlpha2 * sinPhi2)) * bb2 * cc2 +
            ((sinBeta2 * cosAlpha2 * sinPhi2 + cosPhi2 * sinAlpha2) * 
            (sinBeta2 * cosAlpha2 * sinPhi2 + cosPhi2 * sinAlpha2)) * aa2 * cc2 +
            (cosBeta2 * cosBeta2 * cosAlpha2 * cosAlpha2) * aa2 * bb2;
        //xy (must divide by 2f)
        float D2 = ((sinBeta2 * sinAlpha2 * cosPhi2 * cosBeta2 * cosPhi2 + 
            sinPhi2 * cosAlpha2 * cosBeta2 * cosPhi2) * bb2 * cc2 +
            (sinPhi2 * sinBeta2 * sinAlpha2 * cosBeta2 * sinPhi2 - cosAlpha2 * cosPhi2 * cosBeta2 * sinPhi2) * aa2 * cc2 +
            (-sinAlpha2 * cosBeta2 * sinBeta2) * aa2 * bb2);
        //yz (must divide by 2f)
        float E2 = (((-sinBeta2 * cosAlpha2 * cosPhi2 + sinAlpha2 * sinPhi2) * 
            (sinBeta2 * sinAlpha2 * cosPhi2 + sinPhi2 * cosAlpha2)) * bb2 * cc2 +
            ((sinBeta2 * cosAlpha2 * sinPhi2 + cosPhi2 * sinAlpha2) * 
            (-sinPhi2 * sinBeta2 * sinAlpha2 + cosAlpha2 * cosPhi2)) * aa2 * cc2 +
            (-cosBeta2 * cosAlpha2 * cosBeta2 * sinAlpha2) * aa2 * bb2);
        //zx (must divide by 2f)
        float F2 = ((-sinBeta2 * cosAlpha2 * cosPhi2 * cosBeta2 * cosPhi2 + 
            sinAlpha2 * sinPhi2 * cosBeta2 * cosPhi2) * bb2 * cc2 +
            (-sinBeta2 + cosAlpha2 * sinPhi2 * cosBeta2 * sinPhi2 - 
            cosPhi2 * sinAlpha2 * cosBeta2 * sinPhi2) * aa2 * cc2 +
            (cosBeta2 * cosAlpha2 * sinBeta2) * aa2 * bb2);
        //x (must divide by 2f)
        float G2 = - xEllipsoid2 * A2 - yEllipsoid2 * D2 - zEllipsoid2 * F2;
        //y (must divide by 2f)
        float H2 = - xEllipsoid2 * D2 - yEllipsoid2 * B2 - zEllipsoid2 * E2;
        //z (must divide by 2f)
        float I2 = - xEllipsoid2 * F2 - yEllipsoid2 * E2 - zEllipsoid2 * C2;
        //independent
        float J2 = (xEllipsoid2 * xEllipsoid2 * A2) + (xEllipsoid2 * yEllipsoid2 * 2f * D2) + (xEllipsoid2 * zEllipsoid2 * 2f * F2) + (yEllipsoid2 * yEllipsoid2 * B2) +
            (yEllipsoid2 * zEllipsoid2 * 2f * E2) + (zEllipsoid2 * zEllipsoid2 * C2) - (aa2 * bb2 * cc2);

        matrixB[0, 0] = A2;
        matrixB[0, 1] = D2;
        matrixB[0, 2] = F2;
        matrixB[0, 3] = G2;
        matrixB[1, 0] = D2;
        matrixB[1, 1] = B2;
        matrixB[1, 2] = E2;
        matrixB[1, 3] = H2;
        matrixB[2, 0] = F2;
        matrixB[2, 1] = E2;
        matrixB[2, 2] = C2;
        matrixB[2, 3] = I2;
        matrixB[3, 0] = G2;
        matrixB[3, 1] = H2;
        matrixB[3, 2] = I2;
        matrixB[3, 3] = J2;

        //Debug.Log(A2 + " " + B2 + " " + C2 + " " + D2 + " " + E2 + " " + F2 + " " + G2 + " " + H2 + " " + I2 + " " + J2);

        return calcCharacteristicPolynomial(matrixA, matrixB, A1, B1, C1, D1, E1, F1, G1, H1, I1, J1, A2, B2, C2, D2, E2, F2, G2, H2, I2, J2);
    }

    private static float[] calcCharacteristicPolynomial(float[,] matrixA, float[,] matrixB, float A1, float B1, float C1, float D1, float E1, float F1, float G1, float H1, float I1, float J1,
        float A2, float B2, float C2, float D2, float E2, float F2, float G2, float H2, float I2, float J2)
    {
        float a4 = Utility.calDetMatrix4x4(matrixA);
        float a3 = (1 / 3f) * ((A1 * B1 * C1 * J2 + A1 * B1 * C2 * J1 + A1 * B2 * C1 * J1 + A2 * B1 * C1 * J1) +
            (A1 * E1 * I1 * H2 + A1 * E1 * I2 * H1 + A1 * E2 * I1 * H1 + A2 * E1 * I1 * H1) +
            (A1 * H1 * E1 * I2 + A1 * H1 * E2 * I1 + A1 * H2 * E1 * I1 + A2 * H1 * E1 * I1) -
            (A1 * H1 * C1 * H2 + A1 * H1 * C2 * H1 + A1 * H2 * C1 * H1 + A2 * H1 * C1 * H1) -
            (A1 * E1 * E1 * J2 + A1 * E1 * E2 * J1 + A1 * E2 * E1 * J1 + A2 * E1 * E1 * J1) -
            (A1 * B1 * I1 * I2 + A1 * B1 * I2 * I1 + A1 * B2 * I1 * I1 + A2 * B1 * I1 * I1) -
            (D1 * D1 * C1 * J2 + D1 * D1 * C2 * J1 + D1 * D2 * C1 * J1 + D2 * D1 * C1 * J1) -
            (F1 * D1 * I1 * H2 + F1 * D1 * I2 * H1 + F1 * D2 * I1 * H1 + F2 * D1 * I1 * H1) -
            (G1 * D1 * E1 * I2 + G1 * D1 * E2 * I1 + G1 * D2 * E1 * I1 + G2 * D1 * E1 * I1) +
            (G1 * D1 * C1 * H2 + G1 * D1 * C2 * H1 + G1 * D2 * C1 * H1 + G2 * D1 * C1 * H1) +
            (F1 * D1 * E1 * J2 + F1 * D1 * E2 * J1 + F1 * D2 * E1 * J1 + F2 * D1 * E1 * J1) +
            (D1 * D1 * I1 * I2 + D1 * D1 * I2 * I1 + D1 * D2 * I1 * I1 + D2 * D1 * I1 * I1) +
            (D1 * E1 * F1 * J2 + D1 * E1 * F2 * J1 + D1 * E2 * F1 * J1 + D2 * E1 * F1 * J1) +
            (F1 * H1 * F1 * H2 + F1 * H1 * F2 * H1 + F1 * H2 * F1 * H1 + F2 * H1 * F1 * H1) +
            (G1 * B1 * F1 * I2 + G1 * B1 * F2 * I1 + G1 * B2 * F1 * I1 + G2 * B1 * F1 * I1) -
            (G1 * E1 * F1 * H2 + G1 * E1 * F2 * H1 + G1 * E2 * F1 * H1 + G2 * E1 * F1 * H1) -
            (F1 * B1 * F1 * J2 + F1 * B1 * F2 * J1 + F1 * B2 * F1 * J1 + F2 * B1 * F1 * J1) -
            (D1 * H1 * F1 * I2 + D1 * H1 * F2 * I1 + D1 * H2 * F1 * I1 + D2 * H1 * F1 * I1) -
            (D1 * E1 * I1 * G2 + D1 * E1 * I2 * G1 + D1 * E2 * I1 * G1 + D2 * E1 * I1 * G1) -
            (F1 * H1 * E1 * G2 + F1 * H1 * E2 * G1 + F1 * H2 * E1 * G1 + F2 * H1 * E1 * G1) -
            (G1 * B1 * C1 * G2 + G1 * B1 * C2 * G1 + G1 * B2 * C1 * G1 + G2 * B1 * C1 * G1) +
            (G1 * E1 * E1 * G2 + G1 * E1 * E2 * G1 + G1 * E2 * E1 * G1 + G2 * E1 * E1 * G1) +
            (F1 * B1 * I1 * G2 + F1 * B1 * I2 * G1 + F1 * B2 * I1 * G1 + F2 * B1 * I1 * G1) +
            (D1 * H1 * C1 * G2 + D1 * H1 * C2 * G1 + D1 * H2 * C1 * G1 + D2 * H1 * C1 * G1));
        float a2 = (1 / 3f) * ((A1 * B1 * C2 * J2 + A1 * B2 * C1 * J2 + A1 * B2 * C2 * J1 + A2 * B1 * C1 * J2 + A2 * B1 * C2 * J1 + A2 * B2 * C1 * J1) +
            (A1 * E1 * I2 * H2 + A1 * E2 * I1 * H2 + A1 * E2 * I2 * H1 + A2 * E1 * I1 * H2 + A2 * E1 * I2 * H1 + A2 * E2 * I1 * H1) +
            (A1 * H1 * E2 * I2 + A1 * H2 * E1 * I2 + A1 * H2 * E2 * I1 + A2 * H1 * E1 * I2 + A2 * H1 * E2 * I1 + A2 * H2 * E1 * I1) -
            (A1 * H1 * C2 * H2 + A1 * H2 * C1 * H2 + A1 * H2 * C2 * H1 + A2 * H1 * C1 * H2 + A2 * H1 * C2 * H1 + A2 * H2 * C1 * H1) -
            (A1 * E1 * E2 * J2 + A1 * E2 * E1 * J2 + A1 * E2 * E2 * J1 + A2 * E1 * E1 * J2 + A2 * E1 * E2 * J1 + A2 * E2 * E1 * J1) -
            (A1 * B1 * I2 * I2 + A1 * B2 * I1 * I2 + A1 * B2 * I2 * I1 + A2 * B1 * I1 * I2 + A2 * B1 * I2 * I1 + A2 * B2 * I1 * I1) -
            (D1 * D1 * C2 * J2 + D1 * D2 * C1 * J2 + D1 * D2 * C2 * J1 + D2 * D1 * C1 * J2 + D2 * D1 * C2 * J1 + D2 * D2 * C1 * J1) -
            (F1 * D1 * I2 * H2 + F1 * D2 * I1 * H2 + F1 * D2 * I2 * H1 + F2 * D1 * I1 * H2 + F2 * D1 * I2 * H1 + F2 * D2 * I1 * H1) -
            (G1 * D1 * E2 * I2 + G1 * D2 * E1 * I2 + G1 * D2 * E2 * I1 + G2 * D1 * E1 * I2 + G2 * D1 * E2 * I1 + G2 * D2 * E1 * I1) +
            (G1 * D1 * C2 * H2 + G1 * D2 * C1 * H2 + G1 * D2 * C2 * H1 + G2 * D1 * C1 * H2 + G2 * D1 * C2 * H1 + G2 * D2 * C1 * H1) +
            (F1 * D1 * E2 * J2 + F1 * D2 * E1 * J2 + F1 * D2 * E2 * J1 + F2 * D1 * E1 * J2 + F2 * D1 * E2 * J1 + F2 * D2 * E1 * J1) +
            (D1 * D1 * I2 * I2 + D1 * D2 * I1 * I2 + D1 * D2 * I2 * I1 + D2 * D1 * I1 * I2 + D2 * D1 * I2 * I1 + D2 * D2 * I1 * I1) +
            (D1 * E1 * F2 * J2 + D1 * E2 * F1 * J2 + D1 * E2 * F2 * J1 + D2 * E1 * F1 * J2 + D2 * E1 * F2 * J1 + D2 * E2 * F1 * J1) +
            (F1 * H1 * F2 * H2 + F1 * H2 * F1 * H2 + F1 * H2 * F2 * H1 + F2 * H1 * F1 * H2 + F2 * H1 * F2 * H1 + F2 * H2 * F1 * H1) +
            (G1 * B1 * F2 * I2 + G1 * B2 * F1 * I2 + G1 * B2 * F2 * I1 + G2 * B1 * F1 * I2 + G2 * B1 * F2 * I1 + G2 * B2 * F1 * I1) -
            (G1 * E1 * F2 * H2 + G1 * E2 * F1 * H2 + G1 * E2 * F2 * H1 + G2 * E1 * F1 * H2 + G2 * E1 * F2 * H1 + G2 * E2 * F1 * H1) -
            (F1 * B1 * F2 * J2 + F1 * B2 * F1 * J2 + F1 * B2 * F2 * J1 + F2 * B1 * F1 * J2 + F2 * B1 * F2 * J1 + F2 * B2 * F1 * J1) -
            (D1 * H1 * F2 * I2 + D1 * H2 * F1 * I2 + D1 * H2 * F2 * I1 + D2 * H1 * F1 * I2 + D2 * H1 * F2 * I1 + D2 * H2 * F1 * I1) -
            (D1 * E1 * I2 * G2 + D1 * E2 * I1 * G2 + D1 * E2 * I2 * G1 + D2 * E1 * I1 * G2 + D2 * E1 * I2 * G1 + D2 * E2 * I1 * G1) -
            (F1 * H1 * E2 * G2 + F1 * H2 * E1 * G2 + F1 * H2 * E2 * G1 + F2 * H1 * E1 * G2 + F2 * H1 * E2 * G1 + F2 * H2 * E1 * G1) -
            (G1 * B1 * C2 * G2 + G1 * B2 * C1 * G2 + G1 * B2 * C2 * G1 + G2 * B1 * C1 * G2 + G2 * B1 * C2 * G1 + G2 * B2 * C1 * G1) +
            (G1 * E1 * E2 * G2 + G1 * E2 * E1 * G2 + G1 * E2 * E2 * G1 + G2 * E1 * E1 * G2 + G2 * E1 * E2 * G1 + G2 * E2 * E1 * G1) +
            (F1 * B1 * I2 * G2 + F1 * B2 * I1 * G2 + F1 * B2 * I2 * G1 + F2 * B1 * I1 * G2 + F2 * B1 * I2 * G1 + F2 * B2 * I1 * G1) +
            (D1 * H1 * C2 * G2 + D1 * H2 * C1 * G2 + D1 * H2 * C2 * G1 + D2 * H1 * C1 * G2 + D2 * H1 * C2 * G1 + D2 * H2 * C1 * G1));
        float a1 = (1 / 3f) * ((A1 * B2 * C2 * J2 + A2 * B1 * C2 * J2 + A2 * B2 * C1 * J2 + A2 * B2 * C2 * J1) +
            (A1 * E2 * I2 * H2 + A2 * E1 * I2 * H2 + A2 * E2 * I1 * H2 + A2 * E2 * I2 * H1) +
            (A1 * H2 * E2 * I2 + A2 * H1 * E2 * I2 + A2 * H2 * E1 * I2 + A2 * H2 * E2 * I1) -
            (A1 * H2 * C2 * H2 + A2 * H1 * C2 * H2 + A2 * H2 * C1 * H2 + A2 * H2 * C2 * H1) -
            (A1 * E2 * E2 * J2 + A2 * E1 * E2 * J2 + A2 * E2 * E1 * J2 + A2 * E2 * E2 * J1) -
            (A1 * B2 * I2 * I2 + A2 * B1 * I2 * I2 + A2 * B2 * I1 * I2 + A2 * B2 * I2 * I1) -
            (D1 * D2 * C2 * J2 + D2 * D1 * C2 * J2 + D2 * D2 * C1 * J2 + D2 * D2 * C2 * J1) -
            (F1 * D2 * I2 * H2 + F2 * D1 * I2 * H2 + F2 * D2 * I1 * H2 + F2 * D2 * I2 * H1) -
            (G1 * D2 * E2 * I2 + G2 * D1 * E2 * I2 + G2 * D2 * E1 * I2 + G2 * D2 * E2 * I1) +
            (G1 * D2 * C2 * H2 + G2 * D1 * C2 * H2 + G2 * D2 * C1 * H2 + G2 * D2 * C2 * H1) +
            (F1 * D2 * E2 * J2 + F2 * D1 * E2 * J2 + F2 * D2 * E1 * J2 + F2 * D2 * E2 * J1) +
            (D1 * D2 * I2 * I2 + D2 * D1 * I2 * I2 + D2 * D2 * I1 * I2 + D2 * D2 * I2 * I1) +
            (D1 * E2 * F2 * J2 + D2 * E1 * F2 * J2 + D2 * E2 * F1 * J2 + D2 * E2 * F2 * J1) +
            (F1 * H2 * F2 * H2 + F2 * H1 * F2 * H2 + F2 * H2 * F1 * H2 + F2 * H2 * F2 * H1) +
            (G1 * B2 * F2 * I2 + G2 * B1 * F2 * I2 + G2 * B2 * F1 * I2 + G2 * B2 * F2 * I1) -
            (G1 * E2 * F2 * H2 + G2 * E1 * F2 * H2 + G2 * E2 * F1 * H2 + G2 * E2 * F2 * H1) -
            (F1 * B2 * F2 * J2 + F2 * B1 * F2 * J2 + F2 * B2 * F1 * J2 + F2 * B2 * F2 * J1) -
            (D1 * H2 * F2 * I2 + D2 * H1 * F2 * I2 + D2 * H2 * F1 * I2 + D2 * H2 * F2 * I1) -
            (D1 * E2 * I2 * G2 + D2 * E1 * I2 * G2 + D2 * E2 * I1 * G2 + D2 * E2 * I2 * G1) -
            (F1 * H2 * E2 * G2 + F2 * H1 * E2 * G2 + F2 * H2 * E1 * G2 + F2 * H2 * E2 * G1) -
            (G1 * B2 * C2 * G2 + G2 * B1 * C2 * G2 + G2 * B2 * C1 * G2 + G2 * B2 * C2 * G1) +
            (G1 * E2 * E2 * G2 + G2 * E1 * E2 * G2 + G2 * E2 * E1 * G2 + G2 * E2 * E2 * G1) +
            (F1 * B2 * I2 * G2 + F2 * B1 * I2 * G2 + F2 * B2 * I1 * G2 + F2 * B2 * I2 * G1) +
            (D1 * H2 * C2 * G2 + D2 * H1 * C2 * G2 + D2 * H2 * C1 * G2 + D2 * H2 * C2 * G1));
        float a0 = Utility.calDetMatrix4x4(matrixB);

        /*
        Debug.Log("--------------------");
        Debug.Log(a4);
        Debug.Log(a3);
        Debug.Log(a2);
        Debug.Log(a1);
        Debug.Log(a0);
        */

        return new float[5] { a4, a3, a2, a1, a0 };
    }

    public static bool Ellipsoid_Ellipsoid_Caravantes(GameObject Ellipsoid1, GameObject Ellipsoid2)
    {
        float[] characteristicPolynomialValues = characteristicPolynomialEllipsoid(Ellipsoid1, Ellipsoid2);

        float a4 = characteristicPolynomialValues[0];
        float a3 = characteristicPolynomialValues[1];
        float a2 = characteristicPolynomialValues[2];
        float a1 = characteristicPolynomialValues[3];
        float a0 = characteristicPolynomialValues[4];

        /*
        float s0 = a4 * (256f * a0 * a0 * a0 * a4 * a4 * a4 - 192f * a0 * a0 * a1 * a3 * a4 * a4 - 128f * a0 * a0 * a2 * a2 * a4 * a4 + 144f * a0 * a0 *
            a2 * a3 * a3 * a4 - 27f * a0 * a0 * a3 * a3 * a3 * a3 + 144f * a0 * a1 * a1 * a2 * a4 * a4 - 6f * a0 * a1 * a1 * a3 * a3 * a4 - 80f * a0 *
            a1 * a2 * a2 * a3 * a4 + 18f * a0 * a1 * a2 * a3 * a3 * a3 + 16f * a0 * a2 * a2 * a2 * a2 * a4 - 4f * a0 * a2 * a2 * a2 * a3 * a3 - 27f * a1
            * a1 * a1 * a1 * a4 * a4 + 18f * a1 * a1 * a1 * a2 * a3 * a4 - 4f * a1 * a1 * a1 * a3 * a3 * a3 - 4f * a1 * a1 * a2 * a2 * a2 * a4 + a1 * a1
            * a2 * a2 * a3 * a3);
        float s1 = 2f * a4 * (16f * a0 * a2 * a4 * a4 - 6f * a0 * a3 * a3 * a4 - 18f * a1 * a1 * a4 * a4 + 14f * a1 * a2 * a3 * a4 - 3f * a1 * a3 * a3 *
            a3 - 4f * a2 * a2 * a2 * a4 + a2 * a2 * a3 * a3);
        float s10 = -a4 * (48f * a0 * a1 * a4 * a4 - 32f * a0 * a2 * a3 * a4 + 9f * a0 * a3 * a3 * a3 - 3f * a1 * a1 * a3 * a4 + 4f * a1 * a2 * a2 * a4
            - a1 * a2 * a3 * a3);

        Debug.Log("s0: " + s0);
        Debug.Log("s1: " + s1);
        Debug.Log("s10: " + s10);
        */

        //Debug.Log(a4 + " " + a3 + " " + a2 + " " + a1 + " " + a0);

        if (DescartesLawOfSigns.DescartesLawOfSignsFourthDegreePolynomial(a4, a3, a2, a1, a0))
        {
            /*
            if ((s0 == 0f & s1 < 0f & s10 < 0f) || (s0 < 0f))
            {
                return false;
            } else
            {
                return true;
            }
            */
            return false;
        } else
        {
            return true;
        }
    }

    public static float[] characteristicPolynomialEllipsoidEllipticParaboloid(GameObject Ellipsoid1, GameObject EllipticParaboloid)
    {
        GeometryCreator variables1 = Ellipsoid1.GetComponent<GeometryCreator>();
        float xradius1 = variables1.xradius; // Ellipsoid semi axis
        float yradius1 = variables1.yradius;
        float zradius1 = variables1.zradius;
        float xEllipsoid1 = Ellipsoid1.transform.position.x; // Ellipsoid position in space
        float yEllipsoid1 = Ellipsoid1.transform.position.y;
        float zEllipsoid1 = Ellipsoid1.transform.position.z;
        float alpha1 = Utility.getRotationX(Ellipsoid1) * Mathf.Deg2Rad;
        float beta1 = Utility.getRotationY(Ellipsoid1) * Mathf.Deg2Rad;
        float phi1 = Utility.getRotationZ(Ellipsoid1) * Mathf.Deg2Rad;
        float sinAlpha1 = Mathf.Sin(alpha1);
        float sinBeta1 = Mathf.Sin(beta1);
        float sinPhi1 = Mathf.Sin(phi1);
        float cosAlpha1 = Mathf.Cos(alpha1);
        float cosBeta1 = Mathf.Cos(beta1);
        float cosPhi1 = Mathf.Cos(phi1);

        GeometryCreator variables2 = EllipticParaboloid.GetComponent<GeometryCreator>();
        float xradius2 = variables2.xradius; // EllipticParaboloid semi axis
        float yradius2 = variables2.yradius;
        float xEllipticParaboloid2 = EllipticParaboloid.transform.position.x; // EllipticParaboloid position in space
        float yEllipticParaboloid2 = EllipticParaboloid.transform.position.y;
        float zEllipticParaboloid2 = EllipticParaboloid.transform.position.z;
        float alpha2 = Utility.getRotationX(EllipticParaboloid) * Mathf.Deg2Rad;
        float beta2 = Utility.getRotationY(EllipticParaboloid) * Mathf.Deg2Rad;
        float phi2 = Utility.getRotationZ(EllipticParaboloid) * Mathf.Deg2Rad;
        float sinAlpha2 = Mathf.Sin(alpha2);
        float sinBeta2 = Mathf.Sin(beta2);
        float sinPhi2 = Mathf.Sin(phi2);
        float cosAlpha2 = Mathf.Cos(alpha2);
        float cosBeta2 = Mathf.Cos(beta2);
        float cosPhi2 = Mathf.Cos(phi2);

        float[,] matrixA = new float[4, 4];
        float[,] matrixB = new float[4, 4];

        float aa1 = xradius1 * xradius1;
        float bb1 = yradius1 * yradius1;
        float cc1 = zradius1 * zradius1;

        float aa2 = xradius2 * xradius2;
        float bb2 = yradius2 * yradius2;

        //Quadric 1
        //x^2 (/aa1 /bb1 /cc1)
        float A1 = (cosBeta1 * cosBeta1 * cosPhi1 * cosPhi1) * bb1 * cc1 +
            (cosBeta1 * cosBeta1 * sinPhi1 * sinPhi1) * aa1 * cc1 +
            (sinBeta1 * sinBeta1) * aa1 * bb1;
        //y^2
        float B1 = Mathf.Pow((sinBeta1 * sinAlpha1 * cosPhi1 + sinPhi1 * cosAlpha1), 2f) * bb1 * cc1 +
            Mathf.Pow(sinPhi1 * sinBeta1 * sinAlpha1 + cosAlpha1 * cosPhi1, 2f) * aa1 * cc1 +
            Mathf.Pow((sinAlpha1 * cosBeta1), 2f) * aa1 * bb1;
        //z^2
        float C1 = ((-sinBeta1 * cosAlpha1 * cosPhi1 + sinAlpha1 * sinPhi1) * (-sinBeta1 * sinAlpha1 * cosPhi1 +
            sinAlpha1 * sinPhi1)) * bb1 * cc1 +
            ((sinBeta1 * cosAlpha1 * sinPhi1 + cosPhi1 * sinAlpha1) * (sinBeta1 * cosAlpha1 * sinPhi1 +
            cosPhi1 * sinAlpha1)) * aa1 * cc1 +
            (cosBeta1 * cosBeta1 * cosAlpha1 * cosAlpha1) * aa1 * bb1;
        //xy (must divide by 2f)
        float D1 = ((sinBeta1 * sinAlpha1 * cosPhi1 * cosBeta1 * cosPhi1 +
            sinPhi1 * cosAlpha1 * cosBeta1 * cosPhi1) * bb1 * cc1 +
            (sinPhi1 * sinBeta1 * sinAlpha1 * cosBeta1 * sinPhi1 - cosAlpha1 * cosPhi1 * cosBeta1 * sinPhi1) * aa1 * cc1 +
            (-sinAlpha1 * cosBeta1 * sinBeta1) * aa1 * bb1);
        //yz (must divide by 2f)
        float E1 = (((-sinBeta1 * cosAlpha1 * cosPhi1 + sinAlpha1 * sinPhi1) * (sinBeta1 * sinAlpha1 * cosPhi1 +
            sinPhi1 * cosAlpha1)) * bb1 * cc1 +
            ((sinBeta1 * cosAlpha1 * sinPhi1 + cosPhi1 * sinAlpha1) * (-sinPhi1 * sinBeta1 * sinAlpha1 +
            cosAlpha1 * cosPhi1)) * aa1 * cc1 +
            (-cosBeta1 * cosAlpha1 * cosBeta1 * sinAlpha1) * aa1 * bb1);
        //zx (must divide by 2f)
        float F1 = ((-sinBeta1 * cosAlpha1 * cosPhi1 * cosBeta1 * cosPhi1 + sinAlpha1 * sinPhi1 * cosBeta1 * cosPhi1) * bb1 * cc1 +
            (-sinBeta1 + cosAlpha1 * sinPhi1 * cosBeta1 * sinPhi1 - cosPhi1 * sinAlpha1 * cosBeta1 * sinPhi1) * aa1 * cc1 +
            (cosBeta1 * cosAlpha1 * sinBeta1) * aa1 * bb1);

        //x (must divide by 2f)
        float G1 = -xEllipsoid1 * A1 - yEllipsoid1 * D1 - zEllipsoid1 * F1;
        //y (must divide by 2f)
        float H1 = -xEllipsoid1 * D1 - yEllipsoid1 * B1 - zEllipsoid1 * E1;
        //z (must divide by 2f)
        float I1 = -xEllipsoid1 * F1 - yEllipsoid1 * E1 - zEllipsoid1 * C1;
        //independent
        float J1 = (xEllipsoid1 * xEllipsoid1 * A1) + (xEllipsoid1 * yEllipsoid1 * 2f * D1) + (xEllipsoid1 * zEllipsoid1 * 2f * F1) + (yEllipsoid1 * yEllipsoid1 * B1) +
            (yEllipsoid1 * zEllipsoid1 * 2f * E1) + (zEllipsoid1 * zEllipsoid1 * C1) - (aa1 * bb1 * cc1);

        matrixA[0, 0] = A1;
        matrixA[0, 1] = D1;
        matrixA[0, 2] = F1;
        matrixA[0, 3] = G1;
        matrixA[1, 0] = D1;
        matrixA[1, 1] = B1;
        matrixA[1, 2] = E1;
        matrixA[1, 3] = H1;
        matrixA[2, 0] = F1;
        matrixA[2, 1] = E1;
        matrixA[2, 2] = C1;
        matrixA[2, 3] = I1;
        matrixA[3, 0] = G1;
        matrixA[3, 1] = H1;
        matrixA[3, 2] = I1;
        matrixA[3, 3] = J1;

        //Quadric 2
        //x^2
        float A2 = (cosBeta2 * cosBeta2 * cosPhi2 * cosPhi2) * bb2 +
            (cosBeta2 * cosBeta2 * sinPhi2 * sinPhi2) * aa2;
        //y^2
        float B2 = (Mathf.Pow((sinBeta2 * sinAlpha2 * cosPhi2 + sinPhi2 * cosAlpha2), 2f)) * bb2 +
            (Mathf.Pow(sinPhi2 * sinBeta2 * sinAlpha2 + cosAlpha2 * cosPhi2, 2f)) * aa2;
        //z^2
        float C2 = 0f;
        //xy (must divide by 2f)
        float D2 = 0f;
        //yz (must divide by 2f)
        float E2 = 0f;
        //zx (must divide by 2f)
        float F2 = 0f;
        //x (must divide by 2f)
        float G2 = - xEllipticParaboloid2 * A2;
        //y (must divide by 2f)
        float H2 = - yEllipticParaboloid2 * B2;
        //z (must divide by 2f)
        float I2 = -1f;
        //independent
        float J2 = (xEllipticParaboloid2 * xEllipticParaboloid2 * A2) + (yEllipticParaboloid2 * yEllipticParaboloid2 * B2);

        matrixB[0, 0] = A2;
        matrixB[0, 1] = D2;
        matrixB[0, 2] = F2;
        matrixB[0, 3] = G2;
        matrixB[1, 0] = D2;
        matrixB[1, 1] = B2;
        matrixB[1, 2] = E2;
        matrixB[1, 3] = H2;
        matrixB[2, 0] = F2;
        matrixB[2, 1] = E2;
        matrixB[2, 2] = C2;
        matrixB[2, 3] = I2;
        matrixB[3, 0] = G2;
        matrixB[3, 1] = H2;
        matrixB[3, 2] = I2;
        matrixB[3, 3] = J2;

        //Debug.Log(A2 + " " + B2 + " " + C2 + " " + D2 + " " + E2 + " " + F2 + " " + G2 + " " + H2 + " " + I2 + " " + J2);

        return calcCharacteristicPolynomial(matrixA, matrixB, A1, B1, C1, D1, E1, F1, G1, H1, I1, J1, A2, B2, C2, D2, E2, F2, G2, H2, I2, J2);
    }

    public static bool Ellipsoid_EllipticParaboloid_Brozos(GameObject object1, GameObject object2)
    {
        GeometryCreator object1Variables = object1.GetComponent<GeometryCreator>();
        GameObject Ellipsoid;
        GameObject EllipticParaboloid;
        if (object1Variables.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid)
        {
            Ellipsoid = object1;
            EllipticParaboloid = object2;
        }
        else
        {
            Ellipsoid = object2;
            EllipticParaboloid = object1;
        }

        float[] characteristicPolynomialValues = characteristicPolynomialEllipsoidEllipticParaboloid(Ellipsoid, EllipticParaboloid);

        float a4 = characteristicPolynomialValues[0];
        float a3 = characteristicPolynomialValues[1];
        float a2 = characteristicPolynomialValues[2];
        float a1 = characteristicPolynomialValues[3];
        float a0 = characteristicPolynomialValues[4];

        float discriminant = (256f * a0 * a0 * a0 * a4 * a4 * a4 - 192f * a0 * a0 * a1 * a3 * a4 * a4 - 128f * a0 * a0 * a2 * a2 * a4 * a4 + 144f * a0 * a0 *
            a2 * a3 * a3 * a4 - 27f * a0 * a0 * a3 * a3 * a3 * a3 + 144f * a0 * a1 * a1 * a2 * a4 * a4 - 6f * a0 * a1 * a1 * a3 * a3 * a4 - 80f * a0 *
            a1 * a2 * a2 * a3 * a4 + 18f * a0 * a1 * a2 * a3 * a3 * a3 + 16f * a0 * a2 * a2 * a2 * a2 * a4 - 4f * a0 * a2 * a2 * a2 * a3 * a3 - 27f * a1
            * a1 * a1 * a1 * a4 * a4 + 18f * a1 * a1 * a1 * a2 * a3 * a4 - 4f * a1 * a1 * a1 * a3 * a3 * a3 - 4f * a1 * a1 * a2 * a2 * a2 * a4 + a1 * a1
            * a2 * a2 * a3 * a3);

        //Debug.Log(discriminant);
        if (discriminant < 0)
        {
            return true;
        }
        return false;
    }

    public static bool Hyperboloid_Plane(GameObject object1, GameObject object2)
    {
        GeometryCreator object1Variables = object1.GetComponent<GeometryCreator>();
        GameObject Hyperboloid;
        GameObject PlanetoIntersect;
        if (object1Variables.Primitive_3D == GeometryCreator._3Dgeo.Plane)
        {
            PlanetoIntersect = object1;
            Hyperboloid = object2;
        } else
        {
            PlanetoIntersect = object2;
            Hyperboloid = object1;
        }

        if (Hyperboloid.GetComponent<GeometryCreator>().Primitive_3D == GeometryCreator._3Dgeo.OneSurfaceHyperboloid)
        {
            return true;
        }

        MeshFilter filter = PlanetoIntersect.GetComponent<MeshFilter>();
        Vector3 planeNormal = Vector3.zero;

        if (filter && filter.mesh.normals.Length > 0)
            planeNormal = filter.transform.TransformDirection(filter.mesh.normals[0]);

        // Get the Plane Center
        Vector3 planeCenter = PlanetoIntersect.transform.position;

        // Extracting coefficients of the plane equation
        float Ax = planeNormal.x;
        float Ay = planeNormal.y;
        float Az = planeNormal.z;
        float Ad = -Ax * planeCenter.x - Ay * planeCenter.y - Az * planeCenter.z;

        /*
        // Output the coefficients
        Debug.Log("Coefficients of the plane equation: ");
        Debug.Log("Ax: " + Ax);
        Debug.Log("Ay: " + Ay);
        Debug.Log("Az: " + Az);
        Debug.Log("Ad: " + Ad);
        */

        GeometryCreator variables1 = Hyperboloid.GetComponent<GeometryCreator>();
        float c = variables1.zradius;

        if (Mathf.Approximately(Vector3.Dot(planeNormal.normalized, Hyperboloid.transform.forward.normalized), 1f) || 
            Mathf.Approximately(Vector3.Dot(planeNormal.normalized, Hyperboloid.transform.forward.normalized), -1f)) {
            //DrawPlane(Hyperboloid.transform.position+c*Hyperboloid.transform.forward, planeNormal);
            //DrawPlane(Hyperboloid.transform.position-c*Hyperboloid.transform.forward, planeNormal);
            if (IsPlaneBetween(Hyperboloid.transform.position + c * Hyperboloid.transform.forward, PlanetoIntersect.transform.position, Hyperboloid.transform.position - c * Hyperboloid.transform.forward))
            {
                //Debug.Log("There is no collision");
                return false;
            } else
            {
                //Debug.Log("There is collision");
                return true;
            }
        } else if (IsPlaneBetween(Hyperboloid.transform.position + c * Hyperboloid.transform.forward, PlanetoIntersect.transform.position, Hyperboloid.transform.position - c * Hyperboloid.transform.forward))
        {
            float angle = Vector3.Angle(planeNormal, Hyperboloid.transform.forward);
            if (angle > 90f)
            {
                angle = 180f - angle;
            }
            float threshold;
            if (Vector3.Distance(Hyperboloid.transform.position, PlanetoIntersect.transform.position) >= 1f) {
                threshold = (Vector3.Distance(Hyperboloid.transform.position, PlanetoIntersect.transform.position) * c);
            } else
            {
                threshold = c;
            }
            //Debug.Log(angle);
            //Debug.Log(threshold);
            if (angle <= 40f / threshold) {
                //Debug.Log("There is no collision");
                return false;
            }
            else
            {
                //Debug.Log("There is collision");
                return true;
            }
        }
        else
        {
            //Debug.Log("There is collision");
            return true;
        }
    }

    
}
