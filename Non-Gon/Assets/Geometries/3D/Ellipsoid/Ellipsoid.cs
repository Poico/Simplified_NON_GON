using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Ellipsoid : Geometry
{
    public override Vector3 point (float theta, float phi, GeometryCreator variables)
    {
        Vector3 point = new Vector3();

        point[0] = Mathf.Sign(Mathf.Cos(theta))*variables.xradius*Mathf.Abs(Mathf.Cos(theta))*Mathf.Abs(Mathf.Cos(phi));
        point[1] = Mathf.Sign(Mathf.Sin(theta))*variables.yradius*Mathf.Abs(Mathf.Sin(theta))*Mathf.Abs(Mathf.Cos(phi));
        point[2] = Mathf.Sign(Mathf.Sin(phi)) * variables.zradius*Mathf.Abs(Mathf.Sin(phi));

        return point;
    }

    public override Vector3 normal (float theta, float phi, GeometryCreator variables)
    {
        Vector3 normal = new Vector3();

        normal[0] = Mathf.Sign(Mathf.Cos(theta))*(1f/variables.xradius)*Mathf.Abs(Mathf.Cos(theta))*Mathf.Abs(Mathf.Cos(phi));
        normal[1] = Mathf.Sign(Mathf.Sin(theta))*(1f/variables.yradius)*Mathf.Abs(Mathf.Sin(theta))*Mathf.Abs(Mathf.Cos(phi));
        normal[2] = Mathf.Sign(Mathf.Sin(phi))*(1f/variables.zradius)*Mathf.Abs(Mathf.Cos(phi));

        return Vector3.Normalize(normal);
    }

    
}
