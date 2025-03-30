using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Menu : MonoBehaviour
{
    public void Back()
    {
        // Load the previous scene
        int currentSceneIndex = UnityEngine.SceneManagement.SceneManager.GetActiveScene().buildIndex;
        if (currentSceneIndex > 0)
        {
            UnityEngine.SceneManagement.SceneManager.LoadScene(0);
        }
    }

    public void HyperboloidPlane()
    {
        UnityEngine.SceneManagement.SceneManager.LoadScene(1);
    }

    public void CylinderCylinder()
    {
        UnityEngine.SceneManagement.SceneManager.LoadScene(2);
    }

    public void EllipseEllipse2D()
    {
        UnityEngine.SceneManagement.SceneManager.LoadScene(3);
    }

    public void EllipsoidEllipsoid()
    {
        UnityEngine.SceneManagement.SceneManager.LoadScene(4);
    }
    
    public void ConvexCircleCircle()
    {
        UnityEngine.SceneManagement.SceneManager.LoadScene(5);
    }
    public void EllipseEllipse()
    {
        UnityEngine.SceneManagement.SceneManager.LoadScene(6);
    }

    public void JumpingBall()
    {
        UnityEngine.SceneManagement.SceneManager.LoadScene(7);
    }

    public void SuperEllipseLine()
    {
        UnityEngine.SceneManagement.SceneManager.LoadScene(8);
    }
}
