using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Ellipse : MonoBehaviour
{
    public static Vector2 point(float xradius, float yradius, float theta)
    {
        float x = Mathf.Cos (theta) * xradius;
        float y = Mathf.Sin (theta) * yradius;

        return new Vector2(x,y);
    }
}
