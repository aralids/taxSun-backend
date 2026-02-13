def test_tsv_load(test_client, test1):
    filepath = test1["filepath"]
    with open(filepath, "rb") as file:
        response = test_client.post('/load_tsv_data', data={"file": file})
    assert response.status_code == 200