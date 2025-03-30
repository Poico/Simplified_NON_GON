using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public abstract class Geometry
{
    public abstract Vector3 point(float theta, float phi, GeometryCreator variables);
    public abstract Vector3 normal(float theta, float phi, GeometryCreator variables);
}
