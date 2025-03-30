using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;

public class GamesScript : MonoBehaviour
{
    public GameObject object1;
    public GameObject object2;
    public GameObject button;
    public GameObject object3;
    public GameObject object4;
    public GameObject object5;
    public GameObject object6;
    public GameObject object7;
    bool collisionStatus;

    // Update is called once per frame
    void Update()
    {
        if (SceneManager.GetActiveScene().name == "Game1" || SceneManager.GetActiveScene().name == "Game3" || SceneManager.GetActiveScene().name == "Game6")
        {
            collisionStatus = ProximityQuery3D.Ellipsoid_Ellipsoid_Caravantes(object1, object2);
            
        } 
        else if(SceneManager.GetActiveScene().name == "Game2")
        {
            collisionStatus = ProximityQuery2D.Ellipse_Ellipse_Caravantes(object1, object2) || ProximityQuery2D.Ellipse_Ellipse_Caravantes(object1, object3) ||
                ProximityQuery2D.Ellipse_Ellipse_Caravantes(object1, object4) || ProximityQuery2D.Ellipse_Ellipse_Caravantes(object1, object5) ||
                ProximityQuery2D.Ellipse_Ellipse_Caravantes(object1, object6) || ProximityQuery2D.Ellipse_Ellipse_Caravantes(object1, object7);
        }
        else if (SceneManager.GetActiveScene().name == "Game4")
        {
            collisionStatus = ProximityQuery3D.Cylinder_Cylinder_Chittawadigi(object1, object2);
        }
        else if (SceneManager.GetActiveScene().name == "Game5")
        {
            collisionStatus = ProximityQuery3D.Ellipsoid_EllipticParaboloid_Brozos(object1, object2);
        }
        if (collisionStatus)
        {
            button.SetActive(true);
        }
        else
        {
            button.SetActive(false);
        }

    }
}
