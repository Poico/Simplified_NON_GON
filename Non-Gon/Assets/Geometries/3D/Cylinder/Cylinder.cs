using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Cylinder : Geometry
{
    public override Vector3 point (float theta, float phi, GeometryCreator variables)
    {
        Vector3 point = new Vector3();

        point[0] = (Mathf.Cos(theta) * variables.xradius);
        point[1] = (Mathf.Sin(theta) * variables.xradius);
        point[2] = variables.yradius * phi / (Mathf.PI-0.3f) + 0.02f * variables.yradius;

        return point;
    }

    public override Vector3 normal (float theta, float phi, GeometryCreator variables)
    {
        Vector3 point = new Vector3();

        point[0] = (Mathf.Cos(theta) * variables.xradius);
        point[1] = (Mathf.Sin(theta) * variables.xradius);
        point[2] = variables.yradius * phi / (Mathf.PI - 0.3f) + 0.02f * variables.yradius;

        Vector3 normal = Vector3.Normalize(point - variables.transform.position);

        return normal;
    }

}
