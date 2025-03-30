using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Superellipse : MonoBehaviour
{
    public static Vector2 point(float xradius, float yradius, float e, float theta)
    {
        float x;
        float y;
        x = xradius * Mathf.Pow(Mathf.Abs(Mathf.Sin(theta)), e) * Mathf.Sign(Mathf.Sin(theta));
        y = yradius * Mathf.Pow(Mathf.Abs(Mathf.Cos(theta)), e) * Mathf.Sign(Mathf.Cos(theta));
        return new Vector2(x,y);
    }
}
