using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;

public class Draggeable : MonoBehaviour
{
    public GameObject draggeable;

    private void OnMouseDrag()
    {
        Vector3 mousePosition = new Vector3(Input.mousePosition.x, Input.mousePosition.y, -Camera.main.transform.position.z + draggeable.transform.position.z);

        Vector3 draggeablePosition = Camera.main.ScreenToWorldPoint(mousePosition);
        draggeable.transform.position = draggeablePosition;
    }
}
