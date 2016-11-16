
def gen_vert_to_face(num_vert, faces):

    vert_to_face = [[] for i in range(num_vert)]  # a list of empty lists
    for face_ID in range(len(faces)):
        face = faces[face_ID]
        for vert_face_ID in face:
            vert_to_face[vert_face_ID].append(face_ID)

    return vert_to_face
