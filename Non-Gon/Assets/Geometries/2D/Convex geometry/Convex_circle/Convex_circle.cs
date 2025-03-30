using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Convex_circle : MonoBehaviour
{
    public float R = 10f;
    public static float f_c (float alpha, float r)
    {
        float res;
        if(alpha > 0f & alpha < 2f/3f*Mathf.PI)
        {
            res = 2f*r + (1f/3f*Mathf.Pow(alpha,3f)*Mathf.Pow(2f/3f*Mathf.PI - alpha, 4f))*r;
        }
        else
        {
            res = 2f*r;
        }
        return res;
    }

    public static float fd_c (float alpha, float r)
    {
        float res;
        if(alpha > 0f & alpha < 2f/3f*Mathf.PI)
        {
            res = (Mathf.Pow(alpha, 2f)*Mathf.Pow(2f/3f*Mathf.PI - alpha, 4f) - 4f/3f*Mathf.Pow(2f/3f*Mathf.PI - alpha, 3f)*Mathf.Pow(alpha,3f))*r;
        }
        else
        {
            res = 0f;
        }
        return res;
    }

    public static Vector2 point(float angle, float radius)
    {
        float f = f_c(angle, radius);
        float fd = fd_c(angle, radius);
        float xpc = fd*radius/Mathf.Sqrt(f*f + fd*fd);
        float zpc = f*radius/Mathf.Sqrt(f*f + fd*fd) - f;
        float r = Mathf.Sqrt(Mathf.Pow(xpc,2f) + Mathf.Pow(zpc,2f));
        float phi = angle - Mathf.Atan(xpc/zpc);
        Vector2 pos = new Vector2(Mathf.Cos(phi), Mathf.Sin(phi));
        return pos*r;
    }
}
