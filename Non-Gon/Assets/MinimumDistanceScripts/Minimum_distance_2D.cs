using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using Complex = System.Numerics.Complex;
using static Utility;



public class Minimum_distance_2D : MonoBehaviour
{  
    public static Vector2[] Point_Ellipse(Vector2 Point, GameObject Ellipse)
    {
        Vector2[] res= new Vector2[2];
        Vector2 T = new Vector2(); //Closest point on the ellipse 

        GeometryCreator variables = Ellipse.GetComponent<GeometryCreator>();
        float a = variables.xradius; // Ellipse semi axis
        float b = variables.yradius;

        //Tranform the query point to a new system of coordinates relative to the ellipse
        Vector2 Point_ = Ellipse.transform.InverseTransformPoint(Point);
        float s1 = Point_[0];
        float s2 = Point_[1];

        T = Point_Ellipse(Point_, a, b);
        T = Ellipse.transform.TransformPoint(T);

        res[0] = Point;
        res[1] = T;
        
        return res;
    }
    public static Vector2 Point_Ellipse(Vector2 Point, float a, float b)
    {
        Debug.Log("Point_Ellipse");
        Vector2 T = new Vector2(); //Closest point on the ellipse 

        //Tranform the query point to a new system of coordinates relative to the ellipse
        float s1 = Point[0];
        float s2 = Point[1];

        float multy = 1f;
        float multx = 1f;

        if(s2 < 0)
        {
            //Method only works if the quey point is in the fist quadrant, if not, change the point to be in the first quadrant and later flip the closest point
            s2 = -s2;
            multy = -1f;
        }

        if(s2 > 0)
        {
            
            if(s1<0)
            {
                //Method only works if the quey point is in the fist quadrant, if not, change the point to be in the first quadrant and later flip the closest point
                s1 = -s1;
                multx = -1f;
            }
            if(s1>0)
            {
                float a2 = a*a;
                float b2 = b*b;
                float s12 = s1*s1;
                float s22 = s2*s2; 
                float z0 = -a2*a2*b2*b2 + a2*s12*b2*b2 + s22*b2*a2*a2;
                float z1 = 2*b2*s12*a2 + 2*a2*b2*s22 - 2*a2*b2*b2 - 2*b2*a2*a2;
                float z2 = a2*s12 + b2*s22 - b2*b2 - a2*a2 - 4*a2*b2;
                float z3 = -2*(b2+a2);
                //solve quartic equation
                float[] roots = Utility.quartic_roots(-1f,z3,z2,z1,z0);
                //of the four roots, the largest one is the wanted solution
                float t = Mathf.NegativeInfinity;
                for (int r = 0; r < 4; r++)
                {   
                    if (t < roots[r] & ! float.IsNaN(roots[r]))
                    {
                        t = roots[r];
                    }
                }
                T[0] = (a*a*s1)/(t + a*a);
                T[1] = (b*b*s2)/(t + b*b);
            }
            else
            {
                T[0] = 0;
                T[1] = b;
            }
        }
        else
        {
            if(s1 < (a*a - b*b)/a)
            {
                T[0] = a*a*s1/(a*a - b*b);
                T[1] = b*Mathf.Sqrt(1 - (T[0]/a)*(T[0]/a));
            }
            else
            {
                T[0] = a;
                T[1] = 0;
            }
        }
        T[0] = T[0]*multx;
        T[1] = T[1]*multy;
        
        return T;
    }
    public static Vector2[] Ellipse_ellipse (GameObject ellipse1, GameObject ellipse2)
    {
        float tol = 0.1f;
        Vector2[] T = new Vector2[2];
        GeometryCreator variables1 = ellipse1.GetComponent<GeometryCreator>();
        float a1 = variables1.xradius; // Ellipse semi axis
        float b1 = variables1.yradius;
        GeometryCreator variables2 = ellipse2.GetComponent<GeometryCreator>();
        float a2 = variables2.xradius; // Ellipse semi axis
        float b2 = variables2.yradius;

        Vector2 p1 = ellipse1.transform.position;
        Vector2 p2 = Point_Ellipse(p1, ellipse2)[1];
        float dist = Vector2.Distance(p1, p2);
        while(true)
        {
            p1 = Point_Ellipse(p2, ellipse1)[1];
            float dist_ = Vector2.Distance(p1, p2);
            if(Mathf.Abs(dist-dist_) < tol)
            {
                break;
            }
            dist = dist_;
            p2 = Point_Ellipse(p1, ellipse2)[1];
            dist_ = Vector2.Distance(p1, p2);
            if(Mathf.Abs(dist-dist_) < tol)
            {
                break;
            }
            dist = dist_;
        }
        T[0] = p1;
        T[1] = p2;
        return T;
    }
    public static float Ellipse_ellipse_distance_of_closest_approach (float angle1, float angle2, float a1, float a2, float b1, float b2)
    {
        if(b1>a1)
        {
            float a1_ = b1;
            b1 = a1;
            a1 = a1_;
            angle1 = Mathf.PI - angle1;
        }
        if(b2>a2)
        {
            float a2_ = b2;
            b2 = a2;
            a2 = a2_;
            angle2 = Mathf.PI - angle2;
        }

        float e1 = Mathf.Sqrt(1f - (b1*b1)/(a1*a1));
        float e2 = Mathf.Sqrt(1f - (b2*b2)/(a2*a2));

        float kd1 = Mathf.Cos(angle1);
        float kd2 = Mathf.Cos(angle2);

        float kk = Mathf.Cos(angle2-angle1);
        float nu = a1/b1 - 1.0f;

        float[,] Ap = new float[2,2];
        Ap[0,0] = b1*b1/(b2*b2)*(1.0f + 0.5f*(1.0f + kk)*(nu*(2.0f + nu) - e2*e2*(1.0f + nu*kk)*(1.0f + nu*kk)));
        Ap[0,1] = b1*b1/(b2*b2)*0.5f*Mathf.Sqrt(1.0f - kk*kk)*(nu*(2.0f + nu) + e2*e2*(1.0f - nu*nu*kk*kk));
        Ap[1,1] = b1*b1/(b2*b2)*(1.0f + 0.5f*(1.0f - kk)*(nu*(2.0f + nu) - e2*e2*(1.0f - nu*kk)*(1.0f - nu*kk)));

        float lambdap = 0.5f*(Ap[0,0] + Ap[1,1]) + Mathf.Sqrt(0.25f*(Ap[0,0] - Ap[1,1])*(Ap[0,0] - Ap[1,1]) + Ap[0,1]*Ap[0,1]);
        float lambdam = 0.5f*(Ap[0,0] + Ap[1,1]) - Mathf.Sqrt(0.25f*(Ap[0,0] - Ap[1,1])*(Ap[0,0] - Ap[1,1]) + Ap[0,1]*Ap[0,1]);

        float b_2 = 1.0f/Mathf.Sqrt(lambdap);
        float a_2 = 1.0f/Mathf.Sqrt(lambdam);

        float cosphi;
        Complex d_;

        if(kk == 1.0f)
        {
            if(Ap[0,0] > Ap[1,1])
            {
                cosphi = b1/a1*kd1/Mathf.Sqrt(1.0f - e1*e1*kd1*kd1);
            }
            else
            {
                cosphi = Mathf.Sqrt(1.0f - kd1*kd1)/Mathf.Sqrt(1.0f - e1*e1*kd1*kd1);
            }
        }
        else
        {
            cosphi = 1.0f/Mathf.Sqrt(2.0f*(Ap[0,1]*Ap[0,1] + (lambdap - Ap[0,0])*(lambdap - Ap[0,0]))*(1.0f - e1*e1*kd1*kd1))*(Ap[0,1]/Mathf.Sqrt(1.0f + kk)*(b1/a1*kd1 + kd2 + (b1/a1 - 1.0f)*kd1*kk) + (lambdap - Ap[0,0])/Mathf.Sqrt(1.0f - kk)*(b1/a1*kd1 - kd2 - (b1/a1 - 1.0f)*kd1*kk));
        }

        float delta = a_2*a_2/(b_2*b_2) - 1.0f;
        if(delta == 0 || cosphi == 0)
        {
            d_ = 1.0f + a_2;
        }
        else
        {
            float tan2 = 1.0f/(cosphi*cosphi) - 1.0f;
            float A = -(1.0f + tan2)/(b_2*b_2);
            float B = -2.0f*(1.0f + tan2 + delta)/b_2;
            float C = -tan2 - (1.0f + delta)*(1.0f + delta) + (1.0f + (1.0f + delta)*tan2)/(b_2*b_2);
            float D = 2.0f*(1.0f + tan2)*(1.0f + delta)/b_2;
            float E = (1.0f + tan2 + delta)*(1.0f + delta);

            float alpha = -3.0f*B*B/(8.0f*A*A) + C/A;
            float beta = B*B*B/(8.0f*A*A*A) - B*C/(2.0f*A*A) + D/A;
            float gamma = -3.0f*B*B*B*B/(256.0f*A*A*A*A) + C*B*B/(16.0f*A*A*A) - B*D/(4.0f*A*A) + E/A;

            Complex q;
            if(beta == 0)
            {
                q = -B/(4.0f*A) + Complex.Sqrt(0.5f*(-alpha + Complex.Sqrt(alpha*alpha - 4.0f*gamma)));
            }
            else
            {
                float P = -alpha*alpha/12.0f - gamma;
                float Q = -alpha*alpha*alpha/108.0f + alpha*gamma/3.0f - beta*beta/8.0f;
                Complex U = Utility.CubeRootc(-Q*0.5f + Complex.Sqrt(Q*Q*0.25f + P*P*P/27.0f));
                Complex y;
                if(U == 0)
                {
                    y = -5.0f*alpha/6.0f - Utility.CubeRoot(Q);
                }
                else
                {
                    y = -5.0f*alpha/6.0f + U - P/(3.0f*U);
                }
                q = -B/(4.0f*A) + 0.5f*(Complex.Sqrt(alpha + 2.0f*y) + Complex.Sqrt(-(3.0f*alpha + 2.0f*y + 2.0f*beta/Complex.Sqrt(alpha + 2.0f*y))));
            }
            d_ = Complex.Sqrt((q*q - 1.0f)/delta*(1.0f + b_2*(1.0f + delta)/q)*(1.0f + b_2*(1.0f + delta)/q) + (1.0f - (q*q - 1.0f)/delta)*(1.0f + b_2/q)*(1.0f + b_2/q));
        }
        
        return (float)d_.Real*b1/Mathf.Sqrt(1.0f - e1*e1*kd1*kd1);
    }
    public static float Ellipse_ellipse_distance_of_closest_approach (GameObject ellipse1, GameObject ellipse2)
    {
        GeometryCreator variables1 = ellipse1.GetComponent<GeometryCreator>();
        float a1 = variables1.xradius; 
        float b1 = variables1.yradius;
        GeometryCreator variables2 = ellipse2.GetComponent<GeometryCreator>();
        float a2 = variables2.xradius; 
        float b2 = variables2.yradius;
        

        Vector2 k1 = ellipse1.transform.TransformDirection(new Vector2(1,0));
        Vector2 k2 = ellipse2.transform.TransformDirection(new Vector2(1,0));
        Vector2 dl = ellipse1.transform.position - ellipse2.transform.position;
        float angle1 = Vector2.Angle(k1, dl);
        float angle2 = Vector2.Angle(k2, dl);

        float dist = Ellipse_ellipse_distance_of_closest_approach(angle1, angle2, a1, a2, b1, b2);

        return dist;
    }
    public static Vector2[] Superellipse_line (GameObject line, GameObject superellipse)
    {
        GeometryCreator variables = superellipse.GetComponent<GeometryCreator>();
        Vector2 n = line.transform.TransformDirection(new Vector2(0,1));
        float a = variables.xradius; 
        float b = variables.yradius;
        float e = variables.e;

        n = superellipse.transform.InverseTransformDirection(n);
        float nx = n[0];
        float ny = n[1];
        float phi = Mathf.Atan((Mathf.Sign(ny)*Mathf.Pow(Mathf.Abs(b*ny),1f/(2f-e)))/(Mathf.Sign(nx)*Mathf.Pow(Mathf.Abs(a*nx),(1f/(2f-e)))));

        float cos = Mathf.Cos(phi);
        float sin = Mathf.Sin(phi);
        Vector2 T = new Vector2(Mathf.Sign(cos)*a*Mathf.Pow(Mathf.Abs(cos), e), Mathf.Sign(sin)*b*Mathf.Pow(Mathf.Abs(sin), e));
        Vector2 Ti = -T;
        T = superellipse.transform.TransformPoint(T);
        Ti = superellipse.transform.TransformPoint(Ti);
        Vector2 T_ = line.transform.InverseTransformPoint(T);
        Vector2 Ti_ = line.transform.InverseTransformPoint(Ti);

        if(Mathf.Abs(Ti_[1])<Mathf.Abs(T_[1])){
            T = Ti;
            T_ = Ti_;
            }

        Vector2 L = line.transform.TransformPoint(new Vector2(T_[0],0));
        Vector2[] p = new Vector2[2];
        p[0] = L;
        p[1] = T;

        return p;
    } 
    
    public static Vector2[] Convex_Circle (GameObject circle, GameObject Convex)
    {
        Vector2 center_convex = Convex.transform.position;
        GeometryCreator variables = circle.GetComponent<GeometryCreator>();
        float R_c = variables.radius;
        Vector2 center_circle = circle.transform.position;
        Vector2 cc = center_convex - center_circle;
        Vector2[] p = new Vector2[2];
        Vector2 y_ = Convex.transform.TransformDirection(Vector2.left);
        float alpha = Mathf.Deg2Rad*Vector2.SignedAngle(y_, cc);
        Vector2 rpc = Convex_circle.point(alpha, R_c);
        Vector2 convex_point = Convex.transform.TransformPoint(rpc);
        Vector2 l = circle.transform.InverseTransformPoint(convex_point);
        l = l.normalized*R_c;
        l = circle.transform.TransformPoint(l);
        p[0] = l;
        p[1] = convex_point;
        return p;
    }
}
