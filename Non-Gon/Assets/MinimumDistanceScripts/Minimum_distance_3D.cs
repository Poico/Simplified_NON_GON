using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using static Minimum_distance_2D;

public class Minimum_distance_3D : MonoBehaviour
{
    public static Vector3[] point_Ellipsoid (Vector3 point, GameObject Ellipsoid)
    {
        Debug.Log("point_Ellipsoid");
        Vector3[] sol = new Vector3[2];
        GeometryCreator variables = Ellipsoid.GetComponent<GeometryCreator>();
        float a = variables.xradius; // Ellipse semi axis
        float b = variables.yradius;
        float c = variables.zradius;

        //Tranform the query point to a new system of coordinates relative to the ellipse
        Vector3 Point_ = Ellipsoid.transform.InverseTransformPoint(point);
        float s1 = Point_[0];
        float s2 = Point_[1];
        float s3 = Point_[2];

        float x = 0;
        float y = 0;
        float z = 0;

        float multx = 1;
        float multy = 1;
        float multz = 1;

        if(s3 < 0)
        {
            //Method only works if the query point is in the fist quadrant, if not, change the point to be in the first quadrant and later flip the closest point
            s3 = -s3;
            multz = -1f;
        }
        if(s2 < 0)
        {
            //Method only works if the query point is in the fist quadrant, if not, change the point to be in the first quadrant and later flip the closest point
            s2 = -s2;
            multy = -1f;
        }
        if(s1 < 0)
        {
            //Method only works if the query point is in the fist quadrant, if not, change the point to be in the first quadrant and later flip the closest point
            s1 = -s1;
            multx = -1f;
        }

        float emin = Mathf.Min(Mathf.Min(a,b),c);

        if(s3 > 0)
        {
            if(s2 > 0)
            {
                if(s1 > 0)
                {
                    float z0 = s1/a;
                    float z1 = s2/b;
                    float z2 = s3/c;
                    float g = z0*z0 + z1*z1 + z2*z2 - 1f;
                    if(g != 0)
                    {
                        float r0 = a*a/(emin*emin);
                        float r1 = b*b/(emin*emin);
                        float r2 = c*c/(emin*emin);
                        float tbar = Utility.GetRoot(r0, r1, r2, z0, z1, z2, g, 200);
                        x = r0*s1/(tbar + r0);
                        y = r1*s2/(tbar + r1);
                        z = r2*s3/(tbar + r2);
                    }
                    else
                    {
                        x = s1;
                        y = s2;
                        z = s3;
                    }
                    
                }
                else
                {
                    x = 0;
                    Vector2 l = Minimum_distance_2D.Point_Ellipse(new Vector2(s2, s3), b, c);
                    y = l[0];
                    z = l[1];
                }
            }
            else
            {
                y = 0;
                if(s1 > 0)
                {
                    Vector2 l = Minimum_distance_2D.Point_Ellipse(new Vector2(s1, s3), a, c);
                    x = l[0];
                    z = l[1];
                }
                else
                {
                    x = 0;
                    z = c;
                }
            }
        }
        else
        {
            float d = a*a - c*c;
            float d1 = b*b - c*c;
            float n = a*s1;
            float n1 = b*s2;
            bool computed = false;

            if(n < d && n1 < d1)
            {
                float xde = n/d;
                float xde1 = n1/d1;
                float xdesqr = xde*xde;
                float xde1sqr = xde1*xde1;
                float discr = 1f - xdesqr - xde1sqr;

                if(discr > 0)
                {
                    x = a*xde;
                    y = b*xde1;
                    z = c*Mathf.Sqrt(discr);
                    computed = true;
                }
            }
            if(!computed)
            {
                z = 0;
                Vector2 l = Minimum_distance_2D.Point_Ellipse(new Vector2(s1,s2), a, b);
                x = l[0];
                y = l[1];
            }
        }
        x = multx*x;
        y = multy*y;
        z = multz*z;
        sol[0] = point;
        Vector3 pop = new Vector3(x,y,z);
        pop = Ellipsoid.transform.TransformPoint(pop);
        sol[1] = pop;
        return sol;
    }
    //public static Vector3[] Ellipsoid_ellipsoid (GameObject ellipsoid1, GameObject ellipsoid2)
    //{
    //    Debug.Log("Ellipsoid_ellipsoid");
    //    Vector3[] sol = new Vector3[2];
    //    Vector3 d = ellipsoid1.transform.InverseTransformDirection(ellipsoid1.transform.position - ellipsoid2.transform.position);
    //    float angle1 = Vector3.Angle(Vector3.ProjectOnPlane(d, new Vector3(0,0,1)), new Vector3(1,0,0));
    //    float angle2 = Vector3.Angle(d, new Vector3(0,0,1));
    //    GeometryCreator var1 = ellipsoid1.GetComponent<GeometryCreator>();
    //    
    //    Vector3 point0 = ellipsoid1.transform.position;
    //    Vector3 point1 = point_Ellipsoid(point0, ellipsoid2)[1];
    //    float dist0 = Vector3.Distance(point0, point1);
    //    float dist1;
    //    int n_iter = 0;
    //    while(true)
    //    {
    //        point0 = point1;
    //        point1 = point_Ellipsoid(point0, ellipsoid2)[1];
    //        dist1 = Vector3.Distance(point0, point1);
    //        
    //        point0 = point1;
    //        dist0 = dist1;
    //        point1 = point_Ellipsoid(point0, ellipsoid1)[1];
    //        dist1 = Vector3.Distance(point0, point1);
    //        if(Mathf.Abs(dist0 - dist1) < 1E-1f || n_iter > 15)
    //        {
    //            break;
    //        }
    //        dist0 = dist1;
    //        n_iter ++;
    //    }
    //    sol[0] = point1;
    //    sol[1] = point0;
    //
    //    return sol;
    //}
    
   //public static Vector3[] Superellipsoid_Plane (GameObject plane, GameObject superellipsoid)
   //{
   //    Debug.Log("Superellipsoid_Plane");
   //    Vector3[] sol = new Vector3[2];
   //    GeometryCreator variables = superellipsoid.GetComponent<GeometryCreator>();
   //    Vector3 center = superellipsoid.transform.position;
   //    center = superellipsoid.transform.InverseTransformPoint(center);
   //    center = plane.transform.TransformPoint(center);
   //    Vector3 n;
   //    if(center[1] > 0)
   //    {
   //        n = plane.transform.TransformDirection(new Vector3(0,1,0));
   //    }
   //    else
   //    {
   //        n = plane.transform.TransformDirection(new Vector3(0,-1,0));
   //    }//
   //    float a = variables.xradius; // Ellipse semi axis
   //    float b = variables.yradius;
   //    float c = variables.zradius;
   //    float e1 = variables.e1;
   //    float e2 = variables.e2;//
   //    n = superellipsoid.transform.InverseTransformDirection(n);
   //    float nx = n[0];
   //    float ny = n[1];
   //    float nz = n[2];
   //    float EPS = Mathf.Pow(10,-6);
   //    float phi1;
   //    float phi2;
   //    if(Mathf.Abs(nx) <= EPS && Mathf.Abs(ny) <= EPS)
   //    {
   //        phi1 = Mathf.PI/2f;
   //        phi2 = Mathf.Sign(nz)*Mathf.PI/2f;
   //    }
   //    else if(Mathf.Abs(nx) <= EPS && Mathf.Abs(nz) <= EPS)
   //    {
   //        phi1 = Mathf.PI + Mathf.Sign(ny)*Mathf.PI/2f;
   //        phi2 = 0;
   //    }
   //    else if(Mathf.Abs(ny) <= EPS && Mathf.Abs(ny) <= EPS)
   //    {
   //        phi1 = Mathf.PI + Mathf.Sign(nx)*Mathf.PI;
   //        phi2 = 0;
   //    }
   //    else
   //    {
   //        float anx = Mathf.Abs(a*nx);
   //        float bny = Mathf.Abs(b*ny);
   //        phi1 = Mathf.Atan2(Mathf.Sign(ny)*Mathf.Pow(bny,1f/(2f - e1)), (Mathf.Sign(nx)*Mathf.Pow(anx,1f/(2f - e1))));
   //        
   //        if(anx > bny)
   //        {
   //            float Cphi = Mathf.Cos(phi1);
   //            phi2 = Mathf.Atan2((Mathf.Sign(nz)*Mathf.Pow(Mathf.Abs(c*nz*Mathf.Sign(Cphi)*Mathf.Pow(Mathf.Abs(Cphi),2f - e1)),1f/(2f-e2))),(Mathf.Sign(nx)*Mathf.Pow(anx,1f/(2f-e2))));
   //        }
   //        else
   //        {
   //            float Sphi = Mathf.Sin(phi1);
   //            phi2 = Mathf.Atan2((Mathf.Sign(nz)*Mathf.Pow(Mathf.Abs(c*nz*Mathf.Sign(Sphi)*Mathf.Pow(Mathf.Abs(Sphi),2f - e1)),1f/(2f-e2))),(Mathf.Sign(ny)*Mathf.Pow(bny,1f/(2f-e2))));
   //        }
   //    }//
   //    Superellipsoid superellipsoid1 = superellipsoid.GetComponent<Superellipsoid>();
   //    
   //    Vector3 point_sellipsoid1 = superellipsoid1.point(phi1, phi2, variables);
   //    Vector3 point_sellipsoid2 = -point_sellipsoid1;//
   //    Vector3 point1 = superellipsoid.transform.TransformPoint(point_sellipsoid1);
   //    Vector3 point2 = superellipsoid.transform.TransformPoint(point_sellipsoid2);//
   //    Vector3 point_plane1 = plane.transform.InverseTransformPoint(point1);
   //    Vector3 point_plane2 = plane.transform.InverseTransformPoint(point2);//
   //    float dist1 = Mathf.Abs(point_plane1[1]);
   //    float dist2 = Mathf.Abs(point_plane2[1]);//
   //    point_plane1 = new Vector3 (point_plane1[0], 0, point_plane1[2]);
   //    point_plane2 = new Vector3 (point_plane2[0], 0, point_plane2[2]);//
   //    if(dist1 < dist2)
   //    {
   //        sol[0] = plane.transform.TransformPoint(point_plane1);
   //        sol[1] = point1;
   //    }
   //    else
   //    {
   //        sol[0] = plane.transform.TransformPoint(point_plane2);
   //        sol[1] = point2;
   //    }
  //
   //    return sol;
   //    
   //}
   
    //public static Vector3[] AlmostConvexGeometry_Plane (GameObject plane, GameObject geometry)
    //{
    //    Debug.Log("AlmostConvexGeometry_Plane");
    //    Parameters_3d parameters = geometry.GetComponent<Parameters_3d>();
    //    Vector3 n = plane.transform.TransformDirection(new Vector3(0,1,0));
    //    n = -geometry.transform.InverseTransformDirection(n);
    //    float theta = Mathf.Deg2Rad*Vector3.SignedAngle(Vector3.ProjectOnPlane(n, new Vector3(0,0,1)), new Vector3(1,0,0), new Vector3(0,0,-1));
    //    float phi = Mathf.Deg2Rad*Vector3.Angle(n, new Vector3(0,0,1));
    //    Vector3 point = parameters.point(theta, phi, geometry.GetComponent<GeometryCreator>());
    //    point = geometry.transform.TransformPoint(point);
    //    Vector3 point_ = plane.transform.InverseTransformPoint(point);
    //    point_ = new Vector3(point_[0], 0, point_[2]);
    //    point_ = plane.transform.TransformPoint(point_);
    //    Vector3[] res = new Vector3[2];
    //    res[0] = point_;
    //    res[1] = point;
    //    return res;
    //}
}
    