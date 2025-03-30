using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;
using UnityEngine.UI;

public class XY : MonoBehaviour
{
    public GameObject ellipse1;
    public GameObject ellipse2;
    GeometryCreator ellipse1GeometryCreator;
    GeometryCreator ellipse2GeometryCreator;
    public Slider _sliderx1;
    public Slider _slidery1;
    public Slider _sliderrotation1;
    public Slider _sliderx2;
    public Slider _slidery2;
    public Slider _sliderrotation2;

    public void Start()
    {
        ellipse1GeometryCreator = ellipse1.GetComponent<GeometryCreator>();
        ellipse2GeometryCreator = ellipse2.GetComponent<GeometryCreator>();

        _sliderx1.onValueChanged.AddListener((_value) =>
        {
            ellipse1GeometryCreator.xradius = _value;
        });
        _slidery1.onValueChanged.AddListener((_value) =>
        {
            ellipse1GeometryCreator.yradius = _value;
        }); 
        _sliderrotation1.onValueChanged.AddListener((_value) =>
        {
            ellipse1.transform.eulerAngles = new Vector3(
                ellipse1.transform.eulerAngles.x,
                ellipse1.transform.eulerAngles.y,
                _value);
        });

        _sliderx2.onValueChanged.AddListener((_value) =>
        {
            ellipse2GeometryCreator.xradius = _value;
        });
        _slidery2.onValueChanged.AddListener((_value) =>
        {
            ellipse2GeometryCreator.yradius = _value;
        });
        _sliderrotation2.onValueChanged.AddListener((_value) =>
        {
            ellipse2.transform.eulerAngles = new Vector3(
                ellipse2.transform.eulerAngles.x,
                ellipse2.transform.eulerAngles.y,
                _value);
        });
    }

    public void Back()
    {
        SceneManager.LoadScene("XYorXYZ");
    }
    
}
