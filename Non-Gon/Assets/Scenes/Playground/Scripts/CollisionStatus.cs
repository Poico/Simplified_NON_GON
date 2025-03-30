using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

public class CollisionStatus : MonoBehaviour
{
    [SerializeField] private TextMeshProUGUI _sliderText;
    public GameObject object1;
    public GameObject object2;
    public Image image;

    // Update is called once per frame
    void Update()
    {
        GeometryCreator variables1 = object1.GetComponent<GeometryCreator>();
        GeometryCreator variables2 = object2.GetComponent<GeometryCreator>();

        bool collisionStatus = false;
        if (variables1.Primitive_2D == GeometryCreator._2Dgeo.Ellipse)
        {
            collisionStatus = ProximityQuery2D.Ellipse_Ellipse_Caravantes(object1, object2);

        }
        if (variables1.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid & variables2.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid)
        {
            collisionStatus = ProximityQuery3D.Ellipsoid_Ellipsoid_Caravantes(object1, object2);

        }
        if ((variables1.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid & variables2.Primitive_3D == GeometryCreator._3Dgeo.EllipticParaboloid) ||
            (variables1.Primitive_3D == GeometryCreator._3Dgeo.EllipticParaboloid & variables2.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid))
        {
            collisionStatus = ProximityQuery3D.Ellipsoid_EllipticParaboloid_Brozos(object2, object1);

        }
        if (variables1.Primitive_3D == GeometryCreator._3Dgeo.Cylinder & variables2.Primitive_3D == GeometryCreator._3Dgeo.Cylinder)
        {
            collisionStatus = ProximityQuery3D.Cylinder_Cylinder_Chittawadigi(object1, object2);

        }
        if ((variables1.Primitive_3D == GeometryCreator._3Dgeo.TwoSurfaceHyperboloid & variables2.Primitive_3D == GeometryCreator._3Dgeo.Plane) ||
            (variables1.Primitive_3D == GeometryCreator._3Dgeo.Plane & variables2.Primitive_3D == GeometryCreator._3Dgeo.TwoSurfaceHyperboloid) ||
            (variables1.Primitive_3D == GeometryCreator._3Dgeo.Plane & variables2.Primitive_3D == GeometryCreator._3Dgeo.OneSurfaceHyperboloid) ||
            (variables1.Primitive_3D == GeometryCreator._3Dgeo.OneSurfaceHyperboloid & variables2.Primitive_3D == GeometryCreator._3Dgeo.Plane))
        {
            collisionStatus = ProximityQuery3D.Hyperboloid_Plane(object1, object2);
        }

        _sliderText.text = collisionStatus.ToString();
        if (collisionStatus)
        {
            image.color = Color.green;
        } else
        {
            image.color = Color.red;
        }
        
    }
}
