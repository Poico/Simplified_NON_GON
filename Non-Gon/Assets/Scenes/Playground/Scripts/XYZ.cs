using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;
using UnityEngine.UI;

public class XYZ : MonoBehaviour
{
    public GameObject object1;
    public GameObject object2;
    public Slider _sliderx1;
    public Slider _slidery1;
    public Slider _sliderz1;
    public Slider _sliderrotation1;
    public Slider _sliderx2;
    public Slider _slidery2;
    public Slider _sliderz2;
    public Slider _sliderrotation2;
    public Slider _sliderrotation3;
    public Slider _sliderrotation4;
    public Slider _sliderrotation5;

    public void Start()
    {
        GeometryCreator object1GeometryCreator = object1.GetComponent<GeometryCreator>();
        GeometryCreator object2GeometryCreator = object2.GetComponent<GeometryCreator>();

        if (object1GeometryCreator.Primitive_3D == GeometryCreator._3Dgeo.Plane)
        {
            _sliderrotation3.onValueChanged.AddListener((_value) =>
            {
                object1.transform.eulerAngles = new Vector3(
                        _value,
                        object1.transform.eulerAngles.y,
                        object1.transform.eulerAngles.z);
            });
            _sliderrotation4.onValueChanged.AddListener((_value) =>
            {
                object1.transform.eulerAngles = new Vector3(
                        object1.transform.eulerAngles.x,
                        _value,
                        object1.transform.eulerAngles.z);
            });
            _sliderrotation5.onValueChanged.AddListener((_value) =>
            {
                object1.transform.eulerAngles = new Vector3(
                        object1.transform.eulerAngles.x,
                        object1.transform.eulerAngles.y,
                        _value);
            });
        }
        else
        {
            _sliderx1.onValueChanged.AddListener((_value) =>
            {
                object1GeometryCreator.xradius = _value;
            });
            _slidery1.onValueChanged.AddListener((_value) =>
            {
                object1GeometryCreator.yradius = _value;
            });

            if (object1GeometryCreator.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid ||
                object1GeometryCreator.Primitive_3D == GeometryCreator._3Dgeo.TwoSurfaceHyperboloid)
            {
                _sliderz1.onValueChanged.AddListener((_value) =>
                {
                    object1GeometryCreator.zradius = _value;
                });
            }

            if (object1GeometryCreator.Primitive_3D != GeometryCreator._3Dgeo.Cylinder)
            {

                _sliderrotation1.onValueChanged.AddListener((_value) =>
                {
                    object1.transform.eulerAngles = new Vector3(
                        object1.transform.eulerAngles.x,
                        object1.transform.eulerAngles.y,
                        _value);
                });
            }
            else
            {
                _sliderrotation1.onValueChanged.AddListener((_value) =>
                {
                    object1.transform.eulerAngles = new Vector3(
                        90f - _value,
                        -90f,
                        -90f);
                });
            }
        }

        if (object2GeometryCreator.Primitive_3D == GeometryCreator._3Dgeo.Plane)
        {
            _sliderrotation3.onValueChanged.AddListener((_value) =>
            {
                object2.transform.eulerAngles = new Vector3(
                        _value,
                        object2.transform.eulerAngles.y,
                        object2.transform.eulerAngles.z);
            });
            _sliderrotation4.onValueChanged.AddListener((_value) =>
            {
                object2.transform.eulerAngles = new Vector3(
                        object2.transform.eulerAngles.x,
                        _value,
                        object2.transform.eulerAngles.z);
            });
            _sliderrotation5.onValueChanged.AddListener((_value) =>
            {
                object2.transform.eulerAngles = new Vector3(
                        object2.transform.eulerAngles.x,
                        object2.transform.eulerAngles.y,
                        _value);
            });
        }
        else
        {
            _sliderx2.onValueChanged.AddListener((_value) =>
            {
                object2GeometryCreator.xradius = _value;
            });
            _slidery2.onValueChanged.AddListener((_value) =>
            {
                object2GeometryCreator.yradius = _value;
            });

            if (object2GeometryCreator.Primitive_3D == GeometryCreator._3Dgeo.Ellipsoid ||
                object2GeometryCreator.Primitive_3D == GeometryCreator._3Dgeo.TwoSurfaceHyperboloid)
            {
                _sliderz2.onValueChanged.AddListener((_value) =>
                {
                    object2GeometryCreator.zradius = _value;
                });
            }

            if (object2GeometryCreator.Primitive_3D == GeometryCreator._3Dgeo.EllipticParaboloid)
            {
                _sliderrotation2.onValueChanged.AddListener((_value) =>
                {
                    object2.transform.eulerAngles = new Vector3(
                        -90f - _value,
                        -90f,
                        -90f);
                });
            }
            else
            {
                if (object2GeometryCreator.Primitive_3D != GeometryCreator._3Dgeo.Cylinder)
                {
                    _sliderrotation2.onValueChanged.AddListener((_value) =>
                    {
                        object2.transform.eulerAngles = new Vector3(
                            object2.transform.eulerAngles.x,
                            object2.transform.eulerAngles.y,
                            _value);
                    });
                }
                else
                {
                    _sliderrotation2.onValueChanged.AddListener((_value) =>
                    {
                        object2.transform.eulerAngles = new Vector3(
                            90f - _value,
                            -90f,
                            -90f);
                    });
                }
            }
        }
    }

    public void EllipsoidEllipsoid()
    {
        SceneManager.LoadScene("EllipsoidEllipsoid");
    }

    public void CylinderCylinder()
    {
        SceneManager.LoadScene("CylinderCylinder");
    }

    public void EllipsoidEllipticParaboloid()
    {
        SceneManager.LoadScene("EllipsoidEllipticParaboloid");
    }

    public void HyperboloidPlane()
    {
        SceneManager.LoadScene("HyperboloidPlane");
    }

    public void Back()
    {
        SceneManager.LoadScene("XYorXYZ");
    }
    
}
