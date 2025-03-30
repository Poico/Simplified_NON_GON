using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Parameters_3d : Geometry
{
    public static float f (float alpha, float beta)
    {
        float sina = Mathf.Sin(alpha);
        float cosa = Mathf.Cos(alpha);
        float sinb = Mathf.Sin(beta);
        float cosb = Mathf.Cos(beta);
        float res = Mathf.Sqrt(100*sina*sina*cosb*cosb + 900*sina*sina*sinb*sinb + 400*cosa*cosa);

        return res;
    }

    public static float fda (float alpha, float beta)
    {
        float sina = Mathf.Sin(alpha);
        float cosa = Mathf.Cos(alpha);
        float sinb = Mathf.Sin(beta);
        float cosb = Mathf.Cos(beta);
        float res = (900*sinb*sinb + 100*cosb*cosb - 400)*cosa*sina/Mathf.Sqrt(900*sinb*sinb*sina*sina + 100*cosb*cosb*sina*sina + 400*cosa*cosa);
        
        return res;
    }

    public static float fdb (float alpha, float beta)
    {
        float sina = Mathf.Sin(alpha);
        float cosa = Mathf.Cos(alpha);
        float sinb = Mathf.Sin(beta);
        float cosb = Mathf.Cos(beta);
        float res = 800*sina*sina*cosb*sinb/Mathf.Sqrt(900*sinb*sinb*sina*sina + 100*cosb*cosb*sina*sina + 400*cosa*cosa);
        return res;
    }

    public static float fdd (float alpha, float beta)
    {
        float sina = Mathf.Sin(alpha);
        float cosa = Mathf.Cos(alpha);
        float sinb = Mathf.Sin(beta);
        float cosb = Mathf.Cos(beta);
        float res = -800*sina*sina*sinb*sinb/Mathf.Sqrt(900f*sinb*sinb*sina*sina + 100f*cosb*cosb*sina*sina + 400f*cosa*cosa) + 800f*sina*sina*cosb*cosb/Mathf.Sqrt(900*sina*sina*sinb*sinb + 100*sina*sina*cosb*cosb + 400*cosa*cosa) - 800*sina*sina*cosb*sinb*(1800*sina*sina*cosb*sinb - 200*sina*sina*cosb*sinb)/(2f*Mathf.Pow(900f*sina*sina*sinb*sinb + 100*sina*sina*cosb*cosb + 400*cosa*cosa, 3f/2f));
        return res;
    }

    public override Vector3 point (float theta, float phi, GeometryCreator variables)
    {
        Vector3 res = new Vector3();
        Vector3 ephi = Vector3.Normalize(new Vector3(Mathf.Cos(phi)*Mathf.Cos(theta), Mathf.Cos(phi)*Mathf.Sin(theta), -Mathf.Sin(phi)));
        Vector3 etheta = Vector3.Normalize(new Vector3(-Mathf.Sin(theta), Mathf.Cos(theta), 0));
        Vector3 en = Vector3.Normalize(new Vector3(Mathf.Sin(phi)*Mathf.Cos(theta), Mathf.Sin(phi)*Mathf.Sin(theta), Mathf.Cos(phi)));

        if(phi>0f && phi<Mathf.PI)
        {
            res = fda(phi, theta)*ephi + fdb(phi,theta)/Mathf.Sin(phi) * etheta + f(phi, theta)*en;
        }
        else
        {
            if(phi == 0 || phi == Mathf.PI)
            {
            res = fda(phi, theta)*ephi + fdd(phi,theta)/Mathf.Cos(phi) * etheta + f(phi, theta)*en;
            }

        }
        return res;
    }
    public override Vector3 normal(float theta, float phi, GeometryCreator variables)
    {
        return Vector3.Normalize(new Vector3(Mathf.Sin(phi)*Mathf.Cos(theta), Mathf.Sin(phi)*Mathf.Sin(theta), Mathf.Cos(phi)));
    }
}
