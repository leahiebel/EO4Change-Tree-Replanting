# Daily use — running the EO4Change pipeline

Short reference for everyday work after the one-time GEE setup is done.

---

## TL;DR

You don't need to re-authenticate every day. **Just open a terminal and run your scripts.** The gcloud credentials are cached on disk and refresh themselves automatically.

---

## Every-day workflow

### 1. Open a PowerShell terminal in VS Code

`Terminal` menu → `New Terminal`. Make sure the prompt starts in this folder:

```
PS C:\Users\leahi\Documents\DTU\EO4Change\EO4Change-Project\EO4Change-Tree-Replanting>
```

### 2. Run whatever script you need

```powershell
# Generate a config + AOI from a seed point
python make_config.py --lon 9.795 --lat 55.866 --name my-site --buffer 2000 --project daring-pier-498809-k0

# Run the full time-series pipeline (uses config_DK.yaml)
python gee.py
```

That's it.

---

## Where your credentials live

After running `gcloud auth login` and `gcloud auth application-default login` once, two files are stored on your machine and the Earth Engine Python library picks them up automatically:

| File | Purpose |
|---|---|
| `C:\Users\leahi\AppData\Roaming\gcloud\credentials.db` | gcloud CLI login |
| `C:\Users\leahi\AppData\Roaming\gcloud\application_default_credentials.json` | What `ee.Initialize()` reads (this is the important one) |

These hold **refresh tokens** that don't expire under normal use — your Python code silently exchanges them for short-lived access tokens behind the scenes.

---

## When *do* you need to re-authenticate?

Only in these specific situations:

| Symptom | Fix |
|---|---|
| `ee.Initialize()` says `Reauthentication is needed` or `invalid_grant` | `gcloud auth application-default login` |
| You changed your Google account password | `gcloud auth login` and `gcloud auth application-default login` |
| You revoked the token at https://myaccount.google.com/permissions | `gcloud auth application-default login` |
| You haven't used gcloud for ~6 months | `gcloud auth application-default login` |
| You're on a different computer | Full first-time setup again |

Both commands open a browser, you sign in, click Allow — same as the first time.

---

## Quick sanity check

If anything feels off, run this — it confirms gcloud and EE both work:

```powershell
python -c "import ee; ee.Initialize(project='daring-pier-498809-k0'); print('OK:', ee.Number(1).getInfo())"
```

Expected output:

```
OK: 1
```

If you get an error instead, re-run:

```powershell
gcloud auth application-default login
```

---

## Useful gcloud commands

```powershell
# Who am I currently logged in as?
gcloud auth list

# What project does ADC use as the billing/quota project?
gcloud auth application-default print-access-token   # prints a token if creds are fresh

# Change the quota project (only needed once after setup)
gcloud auth application-default set-quota-project daring-pier-498809-k0

# Nuclear option — wipe all auth, start over
gcloud auth revoke --all
```

---

## Where things are

| What | Where |
|---|---|
| GEE project id | `daring-pier-498809-k0` |
| GEE Code Editor (browser) | https://code.earthengine.google.com/ |
| GEE quotas / billing | https://console.cloud.google.com/iam-admin/quotas?project=daring-pier-498809-k0 |
| Authorised app permissions | https://myaccount.google.com/permissions |
